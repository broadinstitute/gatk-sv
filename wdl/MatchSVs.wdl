version 1.0

import "Structs.wdl"
import "SVConcordance.wdl" as svc

workflow MatchSVs {
  input {
    File vcf_a
    File vcf_b

    File contig_list
    String prefix

    File? clustering_config
    File? stratification_config
    Array[String]? track_names
    Array[File]? track_intervals

    File reference_dict

    Float? java_mem_fraction
    String gatk_docker
    String sv_base_mini_docker
  }

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

  call svc.SVConcordance as MatchBtoA{
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

  output {
    File match_vcf_atob = MatchAtoB.concordance_vcf
    File match_vcf_atob_index = MatchAtoB.concordance_vcf_index

    File match_vcf_btoa = MatchBtoA.concordance_vcf
    File match_vcf_btoa_index = MatchBtoA.concordance_vcf_index
  }
}
