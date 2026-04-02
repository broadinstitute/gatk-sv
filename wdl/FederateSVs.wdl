version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svc

workflow FederateSVs {
  input {
    File vcf_a
    File vcf_b

    String prefix_a
    String prefix_b
    String prefix

    File? scores

    File contig_list
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    String? match_svs_complex_additional_args

    Float? java_mem_fraction
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_split_cpx_intervals
    RuntimeAttr? runtime_attr_score_matches
  }

  if (!defined(scores)) {
    call svc.SVConcordance as MatchAtoB {
      input:
        eval_vcf=vcf_a,
        truth_vcf=vcf_b,
        run_match_svs=true,
        output_prefix=prefix + ".AtoB",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        contig_list=contig_list,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        sv_base_mini_docker=sv_base_mini_docker,
        gatk_docker=gatk_docker
    }

    call svc.SVConcordance as MatchBtoA {
      input:
        eval_vcf=vcf_b,
        truth_vcf=vcf_a,
        run_match_svs=true,
        output_prefix=prefix + ".BtoA",
        clustering_config=clustering_config,
        stratification_config=stratification_config,
        track_names=track_names,
        track_intervals=track_intervals,
        contig_list=contig_list,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        sv_base_mini_docker=sv_base_mini_docker,
        gatk_docker=gatk_docker
    }

    call SplitComplexIntervals as SplitComplexIntervalsA {
      input:
        vcf=vcf_a,
        prefix=prefix_a + ".cpx_intervals",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_split_cpx_intervals
    }

    call SplitComplexIntervals as SplitComplexIntervalsB {
      input:
        vcf=vcf_b,
        prefix=prefix_b + ".cpx_intervals",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_split_cpx_intervals
    }

    call svc.SVConcordance as MatchAtoBComplex {
      input:
        eval_vcf=SplitComplexIntervalsA.cpx_intervals_vcf,
        truth_vcf=SplitComplexIntervalsB.cpx_intervals_vcf,
        run_match_svs=true,
        output_prefix=prefix + ".AtoB.cpx",
        additional_args="~{match_svs_complex_additional_args} --match-complex-across-subtypes",
        contig_list=contig_list,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        sv_base_mini_docker=sv_base_mini_docker,
        gatk_docker=gatk_docker
    }

    call svc.SVConcordance as MatchBtoAComplex {
      input:
        eval_vcf=SplitComplexIntervalsB.cpx_intervals_vcf,
        truth_vcf=SplitComplexIntervalsA.cpx_intervals_vcf,
        run_match_svs=true,
        output_prefix=prefix + ".BtoA.cpx",
        additional_args="~{match_svs_complex_additional_args} --match-complex-across-subtypes",
        contig_list=contig_list,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        sv_base_mini_docker=sv_base_mini_docker,
        gatk_docker=gatk_docker
    }

    call ScoreMatches {
      input:
        match_vcf_atob=MatchAtoB.concordance_vcf,
        match_vcf_btoa=MatchBtoA.concordance_vcf,
        match_vcf_atob_cpx=MatchAtoBComplex.concordance_vcf,
        match_vcf_btoa_cpx=MatchBtoAComplex.concordance_vcf,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_score_matches
    }
  }

  File scores_ = select_first([scores, ScoreMatches.scores])

  call SVFederate {
    input:
      vcf_a=vcf_a,
      vcf_b=vcf_b,
      scores=scores_,
      prefix=prefix,
      prefix_a=prefix_a,
      prefix_b=prefix_b,
      reference_fasta=reference_fasta,
      reference_fasta_fai=reference_fasta_fai,
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker
  }

  output {
    File federated_vcf = SVFederate.out
    File federated_vcf_index = SVFederate.out_index
  }
}

task SplitComplexIntervals {
  input {
    File vcf
    String prefix
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File cpx_intervals_vcf = "~{prefix}.splitcpx.sorted.vcf.gz"
    File cpx_intervals_vcf_index = "~{prefix}.splitcpx.sorted.vcf.gz.tbi"
  }

  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/split_cpx_intervals.py" script} \
      --vcf ~{vcf} \
      --out ~{prefix}.splitcpx.vcf.gz

    bcftools sort ~{prefix}.splitcpx.vcf.gz -O z -o ~{prefix}.splitcpx.sorted.vcf.gz
    tabix ~{prefix}.splitcpx.sorted.vcf.gz
  >>>

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

task ScoreMatches {
  input {
    File match_vcf_atob
    File match_vcf_btoa
    File match_vcf_atob_cpx
    File match_vcf_btoa_cpx
    String prefix
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size([match_vcf_atob, match_vcf_btoa, match_vcf_atob_cpx, match_vcf_btoa_cpx], "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File scores = "~{prefix}.scores.tsv"
  }

  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/score_sv_matches.py" script} \
      --a-to-b ~{match_vcf_atob} \
      --b-to-a ~{match_vcf_btoa} \
      --a-to-b-cpx ~{match_vcf_atob_cpx} \
      --b-to-a-cpx ~{match_vcf_btoa_cpx} \
      --out ~{prefix}.scores.tsv
  >>>

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

task SVFederate {
  input {
    File vcf_a
    File vcf_b

    String prefix_a
    String prefix_b

    String prefix
    File scores
    String? additional_args

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf_a: {
                 localization_optional: true
               }
    vcf_b:  {
                 localization_optional: true
               }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size([vcf_a, vcf_b], "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVFederate \
      -A ~{vcf_a} \
      -B ~{vcf_b} \
      --prefix-A ~{prefix_a} \
      --prefix-B ~{prefix_b} \
      --sv-pairs ~{scores} \
      -R ~{reference_fasta} \
      ~{additional_args} \
      -O ~{prefix}.vcf.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

