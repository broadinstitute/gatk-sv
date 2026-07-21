version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ExtraScatter {
  input {
    Array[File] vcf_tarballs
    Int initial_shard_number
    Array[String] batches
    Int records_per_shard_extra_join
    String contig
    String prefix

    File ploidy_table
    String variant_prefix
    File reference_fasta
    File reference_fasta_fai
    File reference_dict
    Float? java_mem_fraction
    String? additional_args

    String sv_pipeline_docker
    String gatk_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_extra_scatter
    RuntimeAttr? runtime_attr_extra_join
    RuntimeAttr? runtime_attr_extra_concat
    RuntimeAttr? runtime_attr_join_vcfs
  }

  if (initial_shard_number == 90) {

    Int n_batches = length(batches)
    scatter (i in range(n_batches)) {
      call ScatterVcf {
        input:
          vcf_tarball = vcf_tarballs[i],
          shard_number = initial_shard_number,
          prefix = batches[i],
          records_per_shard = records_per_shard_extra_join,
          contig=contig,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }

    Int n_shards = ScatterVcf.n_shards[0]
    scatter ( i in range(n_shards) ) {
      # Naively join across batches
      call JoinVcfs as JoinSpecial {
        input:
          vcf_tarballs=ScatterVcf.shards_tarball,
          shard_number = i,
          ploidy_table=ploidy_table,
          output_prefix="~{prefix}.~{initial_shard_number}.extra.~{i}",
          contig=contig,
          variant_prefix="~{variant_prefix}.~{initial_shard_number}.extra.~{i}",
          fast_mode=false,
          pesr_sample_overlap=0,
          pesr_interval_overlap=1,
          pesr_breakend_window=0,
          depth_sample_overlap=0,
          depth_interval_overlap=1,
          depth_breakend_window=0,
          mixed_sample_overlap=0,
          mixed_interval_overlap=1,
          mixed_breakend_window=0,
          additional_args=additional_args,
          reference_fasta=reference_fasta,
          reference_fasta_fai=reference_fasta_fai,
          reference_dict=reference_dict,
          java_mem_fraction=java_mem_fraction,
          gatk_docker=gatk_docker,
          runtime_attr_override=runtime_attr_extra_join
      }
    }

    call MiniTasks.ConcatVcfs {
      input:
        vcfs=JoinSpecial.out,
        vcfs_idx=JoinSpecial.out_index,
        naive=true,
        outfile_prefix="~{prefix}.extra.concat",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_extra_concat
    }
  }

  if (initial_shard_number != 90) {
    call JoinVcfs {
      input:
        vcf_tarballs=vcf_tarballs,
        shard_number = initial_shard_number,
        ploidy_table=ploidy_table,
        output_prefix=prefix,
        contig=contig,
        variant_prefix=variant_prefix,
        fast_mode=false,
        pesr_sample_overlap=0,
        pesr_interval_overlap=1,
        pesr_breakend_window=0,
        depth_sample_overlap=0,
        depth_interval_overlap=1,
        depth_breakend_window=0,
        mixed_sample_overlap=0,
        mixed_interval_overlap=1,
        mixed_breakend_window=0,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        additional_args=additional_args,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_join_vcfs
    }
  }

  output {
    File joined_shard = select_first([JoinVcfs.out, ConcatVcfs.concat_vcf])
    File joined_shard_index = select_first([JoinVcfs.out_index, ConcatVcfs.concat_vcf_idx])
  }

}

task JoinVcfs {
    input {
        Array[File] vcf_tarballs
        Int shard_number

        File ploidy_table
        String output_prefix

        String? contig

        Boolean? fast_mode
        Boolean? omit_members
        Boolean? enable_cnv
        Boolean? default_no_call

        String? algorithm
        String? insertion_length_summary_strategy
        String? breakpoint_summary_strategy
        String? alt_allele_summary_strategy

        Float? defrag_padding_fraction
        Float? defrag_sample_overlap

        Float? depth_sample_overlap
        Float? depth_interval_overlap
        Float? depth_size_similarity
        Int? depth_breakend_window
        Float? mixed_sample_overlap
        Float? mixed_interval_overlap
        Float? mixed_size_similarity
        Int? mixed_breakend_window
        Float? pesr_sample_overlap
        Float? pesr_interval_overlap
        Float? pesr_size_similarity
        Int? pesr_breakend_window

        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        Float? java_mem_fraction
        String? additional_args
        String? variant_prefix

        String gatk_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores: 1,
                                   mem_gb: 3.75,
                                   disk_gb: ceil(10 + size(vcf_tarballs, "GB") * 3 + size(reference_fasta, "GB")),
                                   boot_disk_gb: 10,
                                   preemptible_tries: 3,
                                   max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
    }
    command <<<
        set -euxo pipefail

        function getJavaMem() {
            # get JVM memory in MiB by getting total memory from /proc/meminfo
            # and multiplying by java_mem_fraction
            cat /proc/meminfo \
                | awk -v MEM_FIELD="$1" '{
                    f[substr($1, 1, length($1)-1)] = $2
                } END {
                    printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
                }'
        }
        JVM_MAX_MEM=$(getJavaMem MemTotal)
        echo "JVM memory: $JVM_MAX_MEM"

        # extract VCFs matching the shard number provided
        while IFS= read -r tar_path; do
            tar -xf "$tar_path" --wildcards "*.shard.~{shard_number}.vcf.gz"
        done < ~{write_lines(vcf_tarballs)}

        # create the arguments file listing the VCFs extracted
        ls -1 *.shard.~{shard_number}.vcf.gz | sed 's/^/-V /' > arguments.txt

        # create tabix indexes for vcfs
        while read vcf; do
            tabix $vcf
        done < <(ls -1 *.shard.~{shard_number}.vcf.gz)

        mkdir -p tmp

        gatk --java-options "-Xmx${JVM_MAX_MEM}" SVCluster \
            --arguments_file arguments.txt \
            --output ~{output_prefix}.vcf.gz \
            --ploidy-table ~{ploidy_table} \
            --reference ~{reference_fasta} \
            ~{"-L " + contig} \
            ~{true="--fast-mode" false="" fast_mode} \
            ~{true="--enable-cnv" false="" enable_cnv} \
            ~{true="--omit-members" false="" omit_members} \
            ~{true="--default-no-call" false="" default_no_call} \
            ~{"--variant-prefix " + variant_prefix} \
            ~{"--algorithm " + algorithm} \
            ~{"--defrag-padding-fraction " + defrag_padding_fraction} \
            ~{"--defrag-sample-overlap " + defrag_sample_overlap} \
            ~{"--depth-sample-overlap " + depth_sample_overlap} \
            ~{"--depth-interval-overlap " + depth_interval_overlap} \
            ~{"--depth-size-similarity " + depth_size_similarity} \
            ~{"--depth-breakend-window " + depth_breakend_window} \
            ~{"--mixed-sample-overlap " + mixed_sample_overlap} \
            ~{"--mixed-interval-overlap " + mixed_interval_overlap} \
            ~{"--mixed-size-similarity " + mixed_size_similarity} \
            ~{"--mixed-breakend-window " + mixed_breakend_window} \
            ~{"--pesr-sample-overlap " + pesr_sample_overlap} \
            ~{"--pesr-interval-overlap " + pesr_interval_overlap} \
            ~{"--pesr-size-similarity " + pesr_size_similarity} \
            ~{"--pesr-breakend-window " + pesr_breakend_window} \
            ~{"--insertion-length-summary-strategy " + insertion_length_summary_strategy} \
            ~{"--breakpoint-summary-strategy " + breakpoint_summary_strategy} \
            ~{"--alt-allele-summary-strategy " + alt_allele_summary_strategy} \
            --tmp-dir ${PWD}/tmp \
            ~{additional_args}
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: gatk_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task ScatterVcf {
  input {
    File vcf_tarball
    Int shard_number
    String prefix
    Int records_per_shard
    Int? threads = 1
    String? contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = ceil(10+ 3*size([vcf_tarball], "GB"))
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 2,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    # extract VCF matching the shard number provided
    tar -xf ~{vcf_tarball} --wildcards "*.shard.~{shard_number}.vcf.gz"
    vcf=( *.shard.~{shard_number}.vcf.gz )
    tabix $vcf

    # scatter batch vcf by number of records
    bcftools +scatter $vcf -o . -O z -p "~{prefix}".shard. --threads ~{threads} -n ~{records_per_shard} ~{"-r " + contig}

    # get number of shards (will match across all batches)
    ls -1 "~{prefix}".shard.*.vcf.gz | wc -l | tr -d ' ' > n_shards.txt

    # tarball shards
    tar -czvf "~{prefix}".shards.tar.gz "~{prefix}".shard.*.vcf.gz

  >>>
  output {
    File shards_tarball = "~{prefix}.shards.tar.gz"
    Int n_shards = read_int("n_shards.txt")
  }
}
