version 1.0

import "Structs.wdl"
import "Utils.wdl" as util

workflow ManuallyReviewBalancedSVsPerBatch {
  input {
    String svtype

    File cohort_vcf
    File cohort_vcf_index

    File? batch_manta_tloc_vcf
    File batch_pe_file

    String batch
    File batch_samples

    File generate_pe_tabix_py_script # for development

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_subset_samples
    RuntimeAttr? runtime_attr_combine_tlocs
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_generate_script
    RuntimeAttr? runtime_attr_collect_pe
  }

  call util.SubsetVcfBySamplesList {
    input:
      vcf=cohort_vcf,
      vcf_idx=cohort_vcf_index,
      list_of_samples=batch_samples,
      outfile_name="~{batch}.~{svtype}.vcf.gz",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_subset_samples
  }

  call SvtkVcf2bed {
    input:
      vcf=SubsetVcfBySamplesList.vcf_subset,
      prefix="~{batch}.~{svtype}",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_vcf2bed
  }

  if (defined(batch_manta_tloc_vcf)) {
    call CombineTlocs {
      input:
        ctx_bed=SvtkVcf2bed.vcf2bed_out,
        manta_tloc_vcf=select_first([batch_manta_tloc_vcf]),
        manta_tloc_vcf_index="~{select_first([batch_manta_tloc_vcf])}.tbi",
        prefix=batch,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_combine_tlocs
    }
  }

  call GenerateCpxReviewScript {
    input:
      bed=select_first([CombineTlocs.combined_tloc_bed, SvtkVcf2bed.vcf2bed_out]),
      prefix="~{batch}.~{svtype}",
      batch_pe_file=basename(batch_pe_file),
      sv_pipeline_docker=sv_pipeline_docker,
      generate_pe_tabix_py_script=generate_pe_tabix_py_script,
      runtime_attr_override=runtime_attr_generate_script
  }

  call CollectPEMetrics {
    input:
      prefix="~{batch}.~{svtype}",
      PE_collect_script=GenerateCpxReviewScript.pe_evidence_collection_script,
      batch_pe_file=batch_pe_file,
      batch_pe_file_index=batch_pe_file + ".tbi",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_collect_pe
  }

  output {
    File batch_pe_evidence = CollectPEMetrics.evidence
  }
}


task CombineTlocs {
  input {
    String prefix
    File ctx_bed
    File manta_tloc_vcf
    File manta_tloc_vcf_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 4 * size([ctx_bed, manta_tloc_vcf], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    # manta tloc vcf to bed file
    svtk vcf2bed --info ALL ~{manta_tloc_vcf} stdout > manta_tloc.bed

    # extract manta tloc without a first breakpoint within 100bp of a CTX
    cut -f1-4 manta_tloc.bed > bp1.manta_tloc.bed
    zcat ~{ctx_bed} | cut -f1-4 > bp1.ctx.bed
    bedtools window -w 100 -v -a bp1.manta_tloc.bed -b bp1.ctx.bed | cut -f4 > unique_manta_tloc.ids.txt
    rm bp1.manta_tloc.bed
    rm bp1.ctx.bed

    # extract manta tloc without a second breakpoint within 100bp of a CTX
    awk -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="CHR2")c=i; if($i=="END2")e=i;}} NR>1{ if($0 != "") print $c, $e, $e+1, $4 }' manta_tloc.bed > bp2.manta_tloc.bed
    zcat ~{ctx_bed} | awk -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="CHR2")c=i; if($i=="END2")e=i;}} NR>1{ if($0 != "") print $c, $e, $e+1, $4 }' > bp2.ctx.bed
    bedtools window -w 100 -v -a bp2.manta_tloc.bed -b bp2.ctx.bed | cut -f4 >> unique_manta_tloc.ids.txt
    rm bp2.manta_tloc.bed
    rm bp2.ctx.bed

    # extract unique manta tlocs from bed, shuffle order of columns to match CTX, merge with CTX
    python <<CODE
import gzip
keep_cols = "#chrom start end name svtype samples CHR2 CPX_TYPE CPX_INTERVALS END END2 SOURCE STRANDS SVLEN SVTYPE UNRESOLVED_TYPE AN AC AF".split()
keep_ids = set()
with open("unique_manta_tloc.ids.txt", 'r') as ids:
  for line in ids:
    keep_ids.add(line.strip())
with open("~{prefix}.CTX.combined.bed", 'w') as header:
  header.write("\t".join(keep_cols) + "\n")
with open("merged.bed", 'w') as out:
  with open("manta_tloc.bed", 'r') as manta_tloc:
    keep_colnums = []
    name_colnum = 4
    for line in manta_tloc:
      fields = line.strip().split('\t')
      if line.startswith("#"):
        name_colnum = fields.index("name")
        for col in keep_cols:
          keep_colnums.append(fields.index(col))
      else:
        if fields[name_colnum] in keep_ids:
          out.write("\t".join([fields[x] for x in keep_colnums]) + "\n")
  with gzip.open("~{ctx_bed}", 'rt') as ctx:
    keep_colnums = []
    for line in ctx:
      fields = line.strip().split('\t')
      if line.startswith("#"):
        for col in keep_cols:
          keep_colnums.append(fields.index(col))
      else:
        out.write("\t".join([fields[x] for x in keep_colnums]) + "\n")
CODE
    rm manta_tloc.bed

    # sort (keeping header) and bgzip
    sort -Vk1,1 -k2,2n -k3,3n merged.bed >> ~{prefix}.CTX.combined.bed
    rm merged.bed
    bgzip ~{prefix}.CTX.combined.bed

  >>>

  output {
    File combined_tloc_bed = "~{prefix}.CTX.combined.bed.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task GenerateCpxReviewScript {
    input {
        File bed
        String prefix
        String batch_pe_file
        File generate_pe_tabix_py_script
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 5,
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
      set -euo pipefail

      python ~{generate_pe_tabix_py_script} \
        -i ~{bed} \
        -b ~{batch_pe_file} \
        -p ~{prefix}.pe_review.txt \
        -c collect_PE_evidence.~{prefix}.sh \

    >>>

    output {
        File pe_evidence_collection_script = "collect_PE_evidence.~{prefix}.sh"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SvtkVcf2bed {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = "~{prefix}.vcf2bed.bed.gz"

  # simple record-by-record processing, overhead should be O(1), with disk space usage increased because the operation
  # is copying input into new format
  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    svtk vcf2bed --info ALL ~{vcf} stdout \
      | bgzip -c \
      > "~{output_file}"
  >>>

  output {
    File vcf2bed_out = output_file
  }
}


task CollectPEMetrics {
  input {
    String prefix
    File batch_pe_file
    File batch_pe_file_index
    File PE_collect_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String local_pe_path = basename(batch_pe_file)
  String local_pe_index_path = basename(batch_pe_file_index)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: ceil(30.0 + size(batch_pe_file, "GiB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    mv ~{batch_pe_file} ~{local_pe_path}
    mv ~{batch_pe_file_index} ~{local_pe_index_path}

    if [ $(wc -c < "~{PE_collect_script}") -gt 0 ]; then
      bash ~{PE_collect_script}
    fi

    touch ~{prefix}.pe_review.txt
    bgzip ~{prefix}.pe_review.txt

  >>>

  output {
    File evidence = "~{prefix}.pe_review.txt.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
