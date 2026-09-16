version 1.0

## Runs the MoChA "Prepare Data" step (https://github.com/freeseek/mocha#prepare-data)
## across a large cohort of per-sample array VCFs, each carrying GT/BAF/LRR
## (and ALLELE_A/ALLELE_B/GC) per genotype:
##   1. Merge every per-sample VCF into one joint multi-sample VCF.
##   2. Strip it down to the minimal binary VCF MoChA expects.
##   3. Build the exclusion-sites list (low-divergence segdups, excess
##      heterozygosity, high missingness), optionally merged with a
##      user-supplied exclusion VCF.
##
## Outputs {prefix}.unphased.bcf(+csi) and {prefix}.xcl.bcf(+csi), ready for
## the next MoChA step (phasing).
##
## Reference inputs (per MoChA's docs):
##   dup            tabix-indexed BED of segmental duplications with a
##                   Jukes-Cantor divergence column
##   call_rate_table tab-delimited sample_id \t call_rate (used to drop
##                   samples below call_rate_threshold before computing
##                   ExcHet/F_MISSING for the exclusion list)
##   extra_xcl_vcf  optional additional sites to exclude, merged in

workflow MochaPrepareData {
  input {
    Array[File] vcfs
    File dup
    File dup_idx
    File call_rate_table
    File? extra_xcl_vcf
    String prefix = "mocha"
    Float call_rate_threshold = 0.97
    String docker = "staphb/bcftools:1.19"
    Int merge_cpu = 4
    Int merge_mem_gb = 16
  }

  call MergeVcfs {
    input:
      vcfs = vcfs,
      prefix = prefix,
      cpu = merge_cpu,
      mem_gb = merge_mem_gb,
      docker = docker,
  }

  call PrepareUnphased {
    input:
      merged_vcf = MergeVcfs.merged_vcf,
      merged_vcf_idx = MergeVcfs.merged_vcf_idx,
      prefix = prefix,
      docker = docker,
  }

  call BuildExclusionList {
    input:
      unphased_bcf = PrepareUnphased.unphased_bcf,
      unphased_bcf_idx = PrepareUnphased.unphased_bcf_idx,
      dup = dup,
      dup_idx = dup_idx,
      call_rate_table = call_rate_table,
      call_rate_threshold = call_rate_threshold,
      extra_xcl_vcf = extra_xcl_vcf,
      prefix = prefix,
      docker = docker,
  }

  output {
    File unphased_bcf = PrepareUnphased.unphased_bcf
    File unphased_bcf_idx = PrepareUnphased.unphased_bcf_idx
    File xcl_bcf = BuildExclusionList.xcl_bcf
    File xcl_bcf_idx = BuildExclusionList.xcl_bcf_idx
  }
}

task MergeVcfs {
  input {
    Array[File] vcfs
    String prefix
    Int cpu
    Int mem_gb
    String docker
  }

  Int disk_gb = ceil(size(vcfs, "GB") * 3) + 50

  command <<<
    set -euo pipefail

    mkdir -p normalized
    i=0
    : > filelist.txt
    for f in ~{sep=" " vcfs}; do
      i=$((i + 1))
      out="normalized/sample_${i}.bcf"
      bcftools view --no-version -Ob -o "$out" --write-index "$f"
      echo "$out" >> filelist.txt
    done

    bcftools merge --no-version -Ob -m none --write-index \
      -o ~{prefix}.merged.bcf -l filelist.txt
  >>>

  output {
    File merged_vcf = "~{prefix}.merged.bcf"
    File merged_vcf_idx = "~{prefix}.merged.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{mem_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task PrepareUnphased {
  input {
    File merged_vcf
    File merged_vcf_idx
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(merged_vcf, "GB") * 3) + 20

  command <<<
    set -euo pipefail

    bcftools annotate --no-version -o ~{prefix}.unphased.bcf -Ob --write-index \
      -x ID,QUAL,^INFO/ALLELE_A,^INFO/ALLELE_B,^INFO/AC,^INFO/GC,^FMT/GT,^FMT/BAF,^FMT/LRR \
      ~{merged_vcf}
  >>>

  output {
    File unphased_bcf = "~{prefix}.unphased.bcf"
    File unphased_bcf_idx = "~{prefix}.unphased.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task BuildExclusionList {
  input {
    File unphased_bcf
    File unphased_bcf_idx
    File dup
    File dup_idx
    File call_rate_table
    Float call_rate_threshold
    File? extra_xcl_vcf
    String prefix
    String docker
  }

  Int disk_gb = ceil(size(unphased_bcf, "GB") * 3) + 20

  command <<<
    set -euo pipefail

    awk -F"\t" -v thr="~{call_rate_threshold}" '$2<thr {print $1}' ~{call_rate_table} > samples_xcl_list.txt

    echo '##INFO=<ID=JK,Number=1,Type=Float,Description="Jukes Cantor">' | \
      bcftools annotate --no-version -Ou -a ~{dup} -c CHROM,FROM,TO,JK -h /dev/stdin \
        ~{unphased_bcf} | \
      bcftools view --no-version -Ou -S ^samples_xcl_list.txt | \
      bcftools +fill-tags --no-version -Ou -t ^Y,MT,chrY,chrM -- -t ExcHet,F_MISSING | \
      bcftools view --no-version -Ou -G | \
      bcftools annotate --no-version -o intermediate.xcl.bcf -Ob --write-index \
        -i 'FILTER!="." && FILTER!="PASS" || INFO/JK<.02 || INFO/ExcHet<1e-6 || INFO/F_MISSING>1-~{call_rate_threshold}' \
        -x ^INFO/JK,^INFO/ExcHet,^INFO/F_MISSING

    EXTRA_XCL="~{if defined(extra_xcl_vcf) then extra_xcl_vcf else ""}"
    if [[ -n "$EXTRA_XCL" ]]; then
      bcftools merge --no-version -o ~{prefix}.xcl.bcf -Ob -m none --write-index \
        intermediate.xcl.bcf "$EXTRA_XCL"
    else
      mv intermediate.xcl.bcf ~{prefix}.xcl.bcf
      mv intermediate.xcl.bcf.csi ~{prefix}.xcl.bcf.csi
    fi
  >>>

  output {
    File xcl_bcf = "~{prefix}.xcl.bcf"
    File xcl_bcf_idx = "~{prefix}.xcl.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "8 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
