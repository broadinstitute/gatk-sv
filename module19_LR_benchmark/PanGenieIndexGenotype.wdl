version 1.0

import "Structs.wdl"

import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks

workflow PanGenieIndexGenotype {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File panel_biallelic_vcf
        File panel_biallelic_vcf_idx

        String index_prefix

        Array[File]? pangenie_index_chromosome_graphs
        Array[File]? pangenie_index_chromosome_kmers
        File? pangenie_index_unique_kmers_map
        File? pangenie_index_path_segments_fasta

        File reference_fasta
        File reference_fasta_fai

        Array[String] chromosomes

        Array[File] input_cram_list
        Array[File?] input_crai_list
        Array[String] sample_name_list

        File convert_to_biallelic_script

        Boolean subset_reads = true

        String? pangenie_extra_args

        String pangenie_docker
        String sv_base_pipeline_docker

        File? monitoring_script

        RuntimeAttr? runtime_attr_cram_index
        RuntimeAttr? runtime_attr_pangenie_index
        RuntimeAttr? runtime_attr_process_case_reads
        RuntimeAttr? runtime_attr_process_case_reads_wo_subset
        RuntimeAttr? runtime_attr_pangenie_genotype
        RuntimeAttr? runtime_attr_pangenie_convert
        RuntimeAttr? runtime_attr_convert_bubbles_to_biallelic
    }

    if(!defined(pangenie_index_unique_kmers_map)){
        call LongReadGenotypeTasks.IndexPanGenieRefPanel as index_pangenie_ref_panel {
            input:
                panel_vcf_gz = panel_vcf_gz,
                panel_vcf_gz_tbi = panel_vcf_gz_tbi,
                reference_fasta = reference_fasta,
                chromosomes = chromosomes,
                output_prefix = index_prefix,
                docker_image = pangenie_docker,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_pangenie_index
        }
    }


    File Pindex_unique_kmers_map = select_first([index_pangenie_ref_panel.pangenie_index_unique_kmers_map, pangenie_index_unique_kmers_map])
    File Pindex_path_segments_fasta = select_first([index_pangenie_ref_panel.pangenie_index_path_segments_fasta, pangenie_index_path_segments_fasta])
    Array[File] Pindex_chromosome_kmers = select_first([index_pangenie_ref_panel.pangenie_index_chromosome_kmers, pangenie_index_chromosome_kmers])
    Array[File] Pindex_chromosome_graphs = select_first([index_pangenie_ref_panel.pangenie_index_chromosome_graphs, pangenie_index_chromosome_graphs])


    scatter (i in range(length(sample_name_list))){
        if (!defined(input_crai_list[i])) {
            call LongReadGenotypeTasks.IndexPanGenieCaseReads as index_pangenie_case_reads {
                # TODO we require the alignments to subset by chromosome; change to start from raw reads
                input:
                    input_cram = input_cram_list[i],
                    docker_image = sv_base_pipeline_docker,
                    monitoring_script = monitoring_script,
                    runtime_attr_override = runtime_attr_cram_index
            }
        }

        if (subset_reads) {
          call LongReadGenotypeTasks.PreprocessPanGenieCaseReads as preprocess_pangenie_case_reads {
              # TODO we require the alignments to subset by chromosome; change to start from raw reads
              input:
                  input_cram = input_cram_list[i],
                  input_cram_idx = select_first([input_crai_list[i], index_pangenie_case_reads.cram_idx]),
                  reference_fasta = reference_fasta,
                  reference_fasta_fai = reference_fasta_fai,
                  output_prefix = sample_name_list[i],
                  chromosomes = chromosomes,
                  docker_image = sv_base_pipeline_docker,
                  monitoring_script = monitoring_script,
                  runtime_attr_override = runtime_attr_process_case_reads
          }
        }

        if (!subset_reads) {
            call LongReadGenotypeTasks.PreprocessPanGenieCaseReadsWithoutSubsetting as preprocess_pangenie_case_reads_wo_sub {
                # TODO we require the alignments to subset by chromosome; change to start from raw reads
                input:
                    input_cram = input_cram_list[i],
                    input_cram_idx = select_first([input_crai_list[i], index_pangenie_case_reads.cram_idx]),
                    reference_fasta = reference_fasta,
                    reference_fasta_fai = reference_fasta_fai,
                    output_prefix = sample_name_list[i],
                    docker_image = sv_base_pipeline_docker,
                    monitoring_script = monitoring_script,
                    runtime_attr_override = runtime_attr_process_case_reads_wo_subset
            }
        }

        call LongReadGenotypeTasks.PanGenieGenotype as pangenie_genotype {
            input:
                pangenie_index_chromosome_graphs = Pindex_chromosome_graphs,
                pangenie_index_chromosome_kmers = Pindex_chromosome_kmers,
                pangenie_index_unique_kmers_map = Pindex_unique_kmers_map,
                pangenie_index_path_segments_fasta = Pindex_path_segments_fasta,
                index_prefix = index_prefix,
                input_fasta = select_first([preprocess_pangenie_case_reads.preprocessed_fasta, preprocess_pangenie_case_reads_wo_sub.preprocessed_fasta]),
                sample_name = sample_name_list[i],
                output_prefix = sample_name_list[i],
                docker_image = pangenie_docker,
                extra_args = pangenie_extra_args,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_pangenie_genotype
        }

        call LongReadGenotypeTasks.ConvertBubblesToBiallelic as convert_bubbles_to_biallelic{
          input:
            input_vcf = pangenie_genotype.genotyping_vcf_gz,
            input_vcf_idx = pangenie_genotype.genotyping_vcf_gz_tbi,

            panel_biallelic_vcf = panel_biallelic_vcf,
            panel_biallelic_vcf_idx = panel_biallelic_vcf_idx,

            convert_to_biallelic_script = convert_to_biallelic_script,

            docker_image = sv_base_pipeline_docker,
            runtime_attr_override = runtime_attr_convert_bubbles_to_biallelic
        }
    }

    output{
      Array[File] pangenie_genotyped_vcf = pangenie_genotype.genotyping_vcf_gz
      Array[File] pangenie_genotyped_vcf_idx = pangenie_genotype.genotyping_vcf_gz_tbi
      Array[File] pangenie_genotyped_biallelic_vcf = convert_bubbles_to_biallelic.biallelic_vcf
      Array[File] pangenie_genotyped_biallelic_vcf_idx = convert_bubbles_to_biallelic.biallelic_vcf_idx

    }
}




