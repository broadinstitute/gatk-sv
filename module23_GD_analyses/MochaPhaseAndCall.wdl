version 1.0

## Runs MoChA's "Phase genotypes" and "Call chromosomal alterations" steps
## (https://github.com/freeseek/mocha#phase-genotypes,
##  https://github.com/freeseek/mocha#call-chromosomal-alterations) on the
## unphased.bcf/xcl.bcf produced by MochaPrepareData.wdl.
##
## Phase genotypes:
##   1. Restrict to sites not in the exclusion list (biallelic only) and
##      strip to minimal INFO/AC + FORMAT/GT.
##   2. Phase each chromosome with SHAPEIT5's phase_common against a
##      chromosome-matched reference panel and genetic map.
##   3. Concatenate the phased per-chromosome BCFs and re-annotate the
##      phased GT back onto the full unphased.bcf (preserving all
##      original INFO/FORMAT fields, e.g. BAF/LRR) to produce {prefix}.bcf.
##
## Call chromosomal alterations:
##   Run `bcftools +mocha` on {prefix}.bcf, excluding xcl_bcf sites, to
##   produce the calls table, per-sample stats table, UCSC bed track, and
##   annotated {prefix}.as.bcf.
##
## Deviation from the literal doc script: the doc pre-splits the sites
## BCF into per-chromosome files with `bcftools +scatter` before phasing.
## This WDL instead passes phase_common's own --region flag directly
## against the whole-genome sites BCF (its index makes this an equivalent,
## simple region subset) so that per-chromosome outputs stay in the exact
## same order as the `chromosomes` input array, with no dependence on the
## +scatter plugin or glob() ordering.
##
## Uses three docker images:
##   bcftools_docker  plain bcftools, for PrepareSites/ConcatAndImportPhase.
##                 Default: staphb/bcftools:1.19 (same as MochaPrepareData.wdl).
##   phase_docker  bcftools is NOT enough here: needs SHAPEIT5's phase_common
##                 binary on PATH. Default:
##                 quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2 (official
##                 bioconda/biocontainers build; confirmed it has no bcftools,
##                 hence the separate bcftools_docker above). Its phase_common
##                 binary is installed as `SHAPEIT5_phase_common`, not
##                 `phase_common` -- PhaseChromosome looks for either name.
##   mocha_docker  bcftools built with the MoChA `+mocha` plugin. Default:
##                 liangyu1/bcftools-mocha:v0.0.1 (bcftools 1.20, mocha
##                 plugin build 2024-05-05) -- verified its `+mocha` flags
##                 match this WDL's CallMocha command exactly, and it was
##                 smoke-tested end-to-end against a synthetic phased VCF.
##                 There is no current official MoChA docker image; avoid
##                 cwhelan/mocha:v1.0 (bcftools 1.10.2, mocha plugin from
##                 2020) since its `+mocha` uses older flag names (-r/--rules
##                 instead of -g/--genome) that don't match this WDL's
##                 command as written.
##
## Reference inputs (per MoChA's docs):
##   ref_panel_vcfs/ref_panel_vcf_idxs  one phasing reference panel BCF+idx
##                     per entry in `chromosomes`, same order (e.g. 1000G
##                     high-coverage or HGDP+1kGP, GRCh38)
##   genetic_map       combined shapeit-format genetic map file, gzipped,
##                     columns: chrom (1-22,23=X), pos, rate, cM
##   input_stats       optional sample_id/computed_gender/call_rate table;
##                     MoChA estimates these from the VCF if omitted
##   cnp/mhc_reg/kir_reg  optional regions passed through to `+mocha`
##
## NOTE on verification: `bcftools +mocha` (mocha_docker default image) was
## smoke-tested end-to-end against a synthetic 2-sample phased BCF with
## GT/BAF/LRR/ALLELE_A/ALLELE_B, confirming its flags match CallMocha's
## command and that it runs to completion. SHAPEIT5_phase_common (phase_docker
## default image) was confirmed present and runnable (--help), but
## PhaseChromosome's full command (phase against a real reference panel) has
## NOT been executed end-to-end -- double-check its output against a small
## region before a full production run.

workflow MochaPhaseAndCall {
  input {
    File unphased_bcf
    File unphased_bcf_idx
    File xcl_bcf
    File xcl_bcf_idx

    Array[String] chromosomes = [
      "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
      "chr21", "chr22", "chrX"
    ]
    Array[File] ref_panel_vcfs
    Array[File] ref_panel_vcf_idxs
    File genetic_map

    String genome_assembly = "GRCh38"
    File? input_stats
    File? cnp
    String? mhc_reg
    String? kir_reg

    String prefix = "mocha"
    String bcftools_docker = "staphb/bcftools:1.19"
    String phase_docker = "quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2"
    String mocha_docker = "liangyu1/bcftools-mocha:v0.0.1"

    Int phase_cpu = 4
    Int phase_mem_gb = 16
  }

  call PrepareSites {
    input:
      unphased_bcf = unphased_bcf,
      unphased_bcf_idx = unphased_bcf_idx,
      xcl_bcf = xcl_bcf,
      xcl_bcf_idx = xcl_bcf_idx,
      prefix = prefix,
      docker = bcftools_docker,
  }

  scatter (i in range(length(chromosomes))) {
    call PhaseChromosome {
      input:
        sites_bcf = PrepareSites.sites_bcf,
        sites_bcf_idx = PrepareSites.sites_bcf_idx,
        chrom = chromosomes[i],
        ref_panel_vcf = ref_panel_vcfs[i],
        ref_panel_vcf_idx = ref_panel_vcf_idxs[i],
        genetic_map = genetic_map,
        prefix = prefix,
        cpu = phase_cpu,
        mem_gb = phase_mem_gb,
        docker = phase_docker,
    }
  }

  call ConcatAndImportPhase {
    input:
      phased_chrom_bcfs = PhaseChromosome.phased_chrom_bcf,
      unphased_bcf = unphased_bcf,
      unphased_bcf_idx = unphased_bcf_idx,
      prefix = prefix,
      docker = bcftools_docker,
  }

  call CallMocha {
    input:
      phased_bcf = ConcatAndImportPhase.phased_bcf,
      phased_bcf_idx = ConcatAndImportPhase.phased_bcf_idx,
      xcl_bcf = xcl_bcf,
      xcl_bcf_idx = xcl_bcf_idx,
      genome_assembly = genome_assembly,
      input_stats = input_stats,
      cnp = cnp,
      mhc_reg = mhc_reg,
      kir_reg = kir_reg,
      prefix = prefix,
      docker = mocha_docker,
  }

  output {
    File phased_bcf = ConcatAndImportPhase.phased_bcf
    File phased_bcf_idx = ConcatAndImportPhase.phased_bcf_idx
    File as_bcf = CallMocha.as_bcf
    File as_bcf_idx = CallMocha.as_bcf_idx
    File calls_tsv = CallMocha.calls_tsv
    File stats_tsv = CallMocha.stats_tsv
    File ucsc_bed = CallMocha.ucsc_bed
  }
}

task PrepareSites {
  input {
    File unphased_bcf
    File unphased_bcf_idx
    File xcl_bcf
    File xcl_bcf_idx
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(unphased_bcf, "GB") * 3) + 20

  command <<<
    set -euo pipefail

    bcftools isec --no-version -Ou --complement --exclude "N_ALT>1" --write 1 \
      ~{unphased_bcf} ~{xcl_bcf} | \
      bcftools view --no-version -Ou --min-ac 0 --exclude-uncalled | \
      bcftools annotate --no-version -o ~{prefix}.sites.bcf -Ob --write-index \
        --remove ID,QUAL,^INFO/AC,^FMT/GT
  >>>

  output {
    File sites_bcf = "~{prefix}.sites.bcf"
    File sites_bcf_idx = "~{prefix}.sites.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task PhaseChromosome {
  input {
    File sites_bcf
    File sites_bcf_idx
    String chrom
    File ref_panel_vcf
    File ref_panel_vcf_idx
    File genetic_map
    String prefix
    Int cpu
    Int mem_gb
    String docker
  }

  Int disk_gb = ceil((size(sites_bcf, "GB") + size(ref_panel_vcf, "GB")) * 2) + 20

  command <<<
    set -euo pipefail

    # quay.io/biocontainers/shapeit5 installs this as `SHAPEIT5_phase_common`
    # rather than `phase_common`; accept either name.
    PHASE_BIN="$(command -v phase_common || command -v SHAPEIT5_phase_common || true)"
    if [ -z "$PHASE_BIN" ]; then
      echo "ERROR: could not find phase_common (or SHAPEIT5_phase_common) on PATH inside phase_docker." >&2
      exit 1
    fi

    chr_num=$(echo "~{chrom}" | sed 's/^chr//')
    zcat ~{genetic_map} | sed 's/^23/X/' | awk -v chr="$chr_num" '$1==chr {print $2,$3,$4}' > genetic_map.txt

    "$PHASE_BIN" \
      --thread ~{cpu} \
      --input ~{sites_bcf} \
      --reference ~{ref_panel_vcf} \
      --map genetic_map.txt \
      --region ~{chrom} \
      --output ~{prefix}.~{chrom}.pgt.bcf
  >>>

  output {
    File phased_chrom_bcf = "~{prefix}.~{chrom}.pgt.bcf"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{mem_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task ConcatAndImportPhase {
  input {
    Array[File] phased_chrom_bcfs
    File unphased_bcf
    File unphased_bcf_idx
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(phased_chrom_bcfs, "GB") * 2 + size(unphased_bcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    for f in ~{sep=" " phased_chrom_bcfs}; do
      bcftools index --force "$f"
    done

    bcftools concat --no-version -o ~{prefix}.pgt.bcf -Ob --write-index \
      ~{sep=" " phased_chrom_bcfs}

    bcftools annotate --no-version -o ~{prefix}.bcf -Ob --write-index \
      --annotations ~{prefix}.pgt.bcf --columns -FMT/GT \
      ~{unphased_bcf}
  >>>

  output {
    File phased_bcf = "~{prefix}.bcf"
    File phased_bcf_idx = "~{prefix}.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task CallMocha {
  input {
    File phased_bcf
    File phased_bcf_idx
    File xcl_bcf
    File xcl_bcf_idx
    String genome_assembly
    File? input_stats
    File? cnp
    String? mhc_reg
    String? kir_reg
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(phased_bcf, "GB") * 3) + 20

  command <<<
    set -euo pipefail

    bcftools +mocha \
      --genome ~{genome_assembly} \
      ~{"--input-stats " + input_stats} \
      --no-version \
      --output ~{prefix}.as.bcf \
      --output-type b \
      --variants ^~{xcl_bcf} \
      --calls ~{prefix}.calls.tsv \
      --stats ~{prefix}.stats.tsv \
      --ucsc-bed ~{prefix}.ucsc.bed \
      --write-index \
      ~{"--cnp " + cnp} \
      ~{"--mhc " + mhc_reg} \
      ~{"--kir " + kir_reg} \
      ~{phased_bcf}
  >>>

  output {
    File as_bcf = "~{prefix}.as.bcf"
    File as_bcf_idx = "~{prefix}.as.bcf.csi"
    File calls_tsv = "~{prefix}.calls.tsv"
    File stats_tsv = "~{prefix}.stats.tsv"
    File ucsc_bed = "~{prefix}.ucsc.bed"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
