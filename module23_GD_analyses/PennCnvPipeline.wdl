version 1.0

## PennCNV pipeline for per-sample array VCFs containing BAF/LRR in the FORMAT field.
##
## Steps:
##   1. ExtractSignal   - pulls Name/Chr/Position/LRR/BAF out of each sample's VCF into
##                         PennCNV's expected tab-delimited signal file format.
##   2. CompilePFB       - builds a population B-allele-frequency (PFB) file from ALL
##                         samples' signal files (skip this and supply pfb_file directly
##                         if you already have a PFB for your array platform).
##   3. CalcGCModel      - (optional) builds a GC-content wave-correction model from the
##                         PFB file and a GC reference; only runs if gc_reference_file is provided.
##   4. DetectCNV        - runs PennCNV's HMM-based CNV calling per sample.
##   5. CleanAndFilterCNV - merges fragmented adjacent calls and filters by probe count / length.
##
## NOTES / ASSUMPTIONS:
##   - Confirm your VCF's FORMAT field tag names for LRR/BAF actually match "LRR" and "BAF"
##     (check with `bcftools view -h your.vcf.gz | grep FORMAT`); adjust the bcftools query
##     format string in ExtractSignal if your pipeline names them differently.
##   - The VCF ID field is assumed to contain the marker/probe name used by your PFB/manifest.
##   - `bcftools_docker` needs bcftools installed; `penncnv_docker` needs PennCNV (perl scripts
##     + compiled C binaries) installed.
##   - PennCNV's perl scripts (compile_pfb.pl, cal_gc_snp.pl, detect_cnv.pl, clean_cnv.pl,
##     filter_cnv.pl) and the bundled hhall.hmm model ship inside the PennCNV install, but their
##     directory is not guaranteed to be on $PATH in every image (this is what caused the
##     "compile_pfb.pl: command not found" error). Every PennCNV task below runs a small
##     `find_penncnv_dir` snippet at the start of its command block that locates the install
##     directory at runtime and adds it to $PATH, so this works regardless of the exact internal
##     layout of whichever penncnv_docker image you use -- no hardcoded path required.
##   - If you already have a cohort PFB (and optionally a GC model) rather than building one from
##     scratch, set `precomputed_pfb_file` (and `precomputed_gcmodel_file`) and the workflow will
##     use those instead of running CompilePFB / CalcGCModel.

workflow PennCNVPipeline {
  input {
    Array[File] vcfs
    Array[File] vcf_indices
    Array[String] sample_ids

    # Optional: skip PFB/GC-model generation by supplying your own.
    File? precomputed_pfb_file
    File? precomputed_gcmodel_file
    File? gc_reference_file          # only used if precomputed_gcmodel_file is not supplied

    Int numsnp_filter = 10
    String length_filter = "100k"

    String bcftools_docker = "staphb/bcftools:1.19"
    String penncnv_docker = "genomicslab/penncnv:1.0.5"
  }

  scatter (idx in range(length(vcfs))) {
    call ExtractSignal {
      input:
        vcf        = vcfs[idx],
        vcf_idx    = vcf_indices[idx],
        sample_id  = sample_ids[idx],
        docker     = bcftools_docker
    }
  }

  if (!defined(precomputed_pfb_file)) {
    call CompilePFB {
      input:
        signal_files = ExtractSignal.signal_file,
        docker       = penncnv_docker
    }
  }

  File pfb_file_final = select_first([precomputed_pfb_file, CompilePFB.pfb_file])

  if (!defined(precomputed_gcmodel_file) && defined(gc_reference_file)) {
    call CalcGCModel {
      input:
        pfb_file          = pfb_file_final,
        gc_reference_file = select_first([gc_reference_file]),
        docker            = penncnv_docker
    }
  }

  File? gcmodel_file_final = if defined(precomputed_gcmodel_file)
                              then precomputed_gcmodel_file
                              else CalcGCModel.gcmodel_file

  scatter (idx in range(length(vcfs))) {
    call DetectCNV {
      input:
        signal_file  = ExtractSignal.signal_file[idx],
        pfb_file     = pfb_file_final,
        gcmodel_file = gcmodel_file_final,
        sample_id    = sample_ids[idx],
        docker       = penncnv_docker
    }

    call CleanAndFilterCNV {
      input:
        rawcnv_file    = DetectCNV.rawcnv_file,
        signal_file    = ExtractSignal.signal_file[idx],
        sample_id      = sample_ids[idx],
        numsnp_filter  = numsnp_filter,
        length_filter  = length_filter,
        docker         = penncnv_docker
    }
  }

  output {
    File pfb_file                    = pfb_file_final
    File? gcmodel_file                = gcmodel_file_final
    Array[File] signal_files          = ExtractSignal.signal_file
    Array[File] raw_cnv_calls         = DetectCNV.rawcnv_file
    Array[File] penncnv_logs          = DetectCNV.log_file
    Array[File] filtered_cnv_calls    = CleanAndFilterCNV.filtered_cnv_file
  }
}

task ExtractSignal {
  input {
    File vcf
    File vcf_idx
    String sample_id
    String docker
  }

  command <<<
    set -euo pipefail
    echo -e "Name\tChr\tPosition\t~{sample_id}.Log R Ratio\t~{sample_id}.B Allele Freq" > ~{sample_id}.pennsignal.txt
    bcftools query -f '%ID\t%CHROM\t%POS[\t%LRR\t%BAF]\n' ~{vcf} >> ~{sample_id}.pennsignal.txt
  >>>

  output {
    File signal_file = "~{sample_id}.pennsignal.txt"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}

task CompilePFB {
  input {
    Array[File] signal_files
    String docker
    String output_basename = "cohort"
  }

  command <<<
    set -euo pipefail
    # Locate compile_pfb.pl wherever it lives inside this image and add it to PATH.
    PENNCNV_DIR="$(dirname "$(find / -xdev -iname 'compile_pfb.pl' 2>/dev/null | head -n1)")"
    if [ -z "$PENNCNV_DIR" ] || [ "$PENNCNV_DIR" = "." ]; then
      echo "ERROR: could not locate compile_pfb.pl inside the container image." >&2
      exit 1
    fi
    export PATH="$PENNCNV_DIR:$PATH"

    compile_pfb.pl -listfile ~{write_lines(signal_files)} -output ~{output_basename}.pfb
  >>>

  output {
    File pfb_file = "~{output_basename}.pfb"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
  }
}

task CalcGCModel {
  input {
    File pfb_file
    File gc_reference_file
    String docker
    String output_basename = "cohort"
  }

  command <<<
    set -euo pipefail
    PENNCNV_DIR="$(dirname "$(find / -xdev -iname 'cal_gc_snp.pl' 2>/dev/null | head -n1)")"
    if [ -z "$PENNCNV_DIR" ] || [ "$PENNCNV_DIR" = "." ]; then
      echo "ERROR: could not locate cal_gc_snp.pl inside the container image." >&2
      exit 1
    fi
    export PATH="$PENNCNV_DIR:$PATH"

    cal_gc_snp.pl ~{gc_reference_file} ~{pfb_file} -output ~{output_basename}.gcmodel
  >>>

  output {
    File gcmodel_file = "~{output_basename}.gcmodel"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}

task DetectCNV {
  input {
    File signal_file
    File pfb_file
    File? gcmodel_file
    String sample_id
    String docker
  }

  command <<<
    set -euo pipefail
    PENNCNV_DIR="$(dirname "$(find / -xdev -iname 'detect_cnv.pl' 2>/dev/null | head -n1)")"
    if [ -z "$PENNCNV_DIR" ] || [ "$PENNCNV_DIR" = "." ]; then
      echo "ERROR: could not locate detect_cnv.pl inside the container image." >&2
      exit 1
    fi
    export PATH="$PENNCNV_DIR:$PATH"

    # hhall.hmm ships inside the same PennCNV install directory tree.
    HMM_PATH="$(find "$PENNCNV_DIR" -iname 'hhall.hmm' 2>/dev/null | head -n1)"
    if [ -z "$HMM_PATH" ]; then
      # fall back to a broader search in case lib/ isn't directly under PENNCNV_DIR
      HMM_PATH="$(find / -xdev -iname 'hhall.hmm' 2>/dev/null | head -n1)"
    fi
    if [ -z "$HMM_PATH" ]; then
      echo "ERROR: could not locate hhall.hmm inside the container image." >&2
      exit 1
    fi

    detect_cnv.pl -test \
      -hmm "$HMM_PATH" \
      -pfb ~{pfb_file} \
      ~{"-gcmodel " + gcmodel_file} \
      -log ~{sample_id}.penncnv.log \
      -out ~{sample_id}.rawcnv \
      ~{signal_file}
  >>>

  output {
    File rawcnv_file = "~{sample_id}.rawcnv"
    File log_file    = "~{sample_id}.penncnv.log"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}

task CleanAndFilterCNV {
  input {
    File rawcnv_file
    File signal_file
    String sample_id
    Int numsnp_filter
    String length_filter
    String docker
  }

  command <<<
    set -euo pipefail
    PENNCNV_DIR="$(dirname "$(find / -xdev -iname 'clean_cnv.pl' 2>/dev/null | head -n1)")"
    if [ -z "$PENNCNV_DIR" ] || [ "$PENNCNV_DIR" = "." ]; then
      echo "ERROR: could not locate clean_cnv.pl inside the container image." >&2
      exit 1
    fi
    export PATH="$PENNCNV_DIR:$PATH"

    clean_cnv.pl combineseg ~{rawcnv_file} -signalfile ~{signal_file} > ~{sample_id}.merged.rawcnv
    filter_cnv.pl ~{sample_id}.merged.rawcnv -numsnp ~{numsnp_filter} -length ~{length_filter} -out ~{sample_id}.filtered.rawcnv
  >>>

  output {
    File filtered_cnv_file = "~{sample_id}.filtered.rawcnv"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}
