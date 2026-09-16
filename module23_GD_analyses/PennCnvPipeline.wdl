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
##   - The HMM model file (hhall.hmm) ships bundled inside the PennCNV install itself, so it is
##     NOT a workflow input -- DetectCNV references it directly at `hmm_path` inside the
##     container. Verify this path matches your chosen penncnv_docker image by running:
##       docker run --rm <penncnv_docker> find / -name "*.hmm" 2>/dev/null
##     and update the `hmm_path` default in the DetectCNV task below if it differs.
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

    # HMM model ships bundled inside the PennCNV install -- not exposed as a workflow input.
    # Verify/update this path for whichever penncnv_docker image you actually use.
    String hmm_path = "/opt/PennCNV/lib/hhall.hmm"
  }

  command <<<
    set -euo pipefail
    detect_cnv.pl -test \
      -hmm ~{hmm_path} \
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
