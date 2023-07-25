# Workflow for running GATK GermlineCNVCaller on multiple case samples using a trained model (obtained from running
# GATK GermlineCNVCaller in the cohort mode). Supports both WGS and WES.
#
# Notes:
#
# - The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
#   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   These intervals will be padded on both sides by the amount specified by padding (default 250)
#   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
#   e.g., for WES).  For WGS, the intervals should simply cover the chromosomes of interest.
#
# - Intervals can be excluded from coverage collection and all downstream steps by using the exclude_intervals
#   argument, which accepts formats compatible with the GATK -XL argument
#   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
#   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
#   be manually filtered from the final callset.
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_case_workflow.wdl -i my_parameters.json
#
#############

version 1.0

import "GermlineCNVTasks.wdl" as CNVTasks

workflow CNVGermlineCaseWorkflow {

    input {
      ##################################
      #### required basic arguments ####
      ##################################
      Array[File] counts
      Array[String] count_entity_ids
      File contig_ploidy_model_tar
      Array[File] gcnv_model_tars
      String gatk_docker
      String linux_docker
      String sv_base_mini_docker

      ##################################
      #### optional basic arguments ####
      ##################################
      File? gatk4_jar_override

      ######################################################################
      #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
      ######################################################################
      Float? ploidy_mapping_error_rate
      Float? ploidy_sample_psi_scale

      ##########################################################
      #### optional arguments for GermlineCNVCallerCaseMode ####
      ##########################################################
      Float? gcnv_p_alt
      Float? gcnv_cnv_coherence_length
      Int? gcnv_max_copy_number

      # optional arguments for germline CNV denoising model
      Float? gcnv_mapping_error_rate
      Float? gcnv_sample_psi_scale
      Float? gcnv_depth_correction_tau
      String? gcnv_copy_number_posterior_expectation_mode
      Int? gcnv_active_class_padding_hybrid_mode

      # optional arguments for Hybrid ADVI
      Float? gcnv_learning_rate
      Float? gcnv_adamax_beta_1
      Float? gcnv_adamax_beta_2
      Int? gcnv_log_emission_samples_per_round
      Float? gcnv_log_emission_sampling_median_rel_error
      Int? gcnv_log_emission_sampling_rounds
      Int? gcnv_max_advi_iter_first_epoch
      Int? gcnv_max_advi_iter_subsequent_epochs
      Int? gcnv_min_training_epochs
      Int? gcnv_max_training_epochs
      Float? gcnv_initial_temperature
      Int? gcnv_num_thermal_advi_iters
      Int? gcnv_convergence_snr_averaging_window
      Float? gcnv_convergence_snr_trigger_threshold
      Int? gcnv_convergence_snr_countdown_window
      Int? gcnv_max_calling_iters
      Float? gcnv_caller_update_convergence_threshold
      Float? gcnv_caller_internal_admixing_rate
      Float? gcnv_caller_external_admixing_rate
      Boolean? gcnv_disable_annealing

      ###################################################
      #### arguments for PostprocessGermlineCNVCalls ####
      ###################################################
      Int ref_copy_number_autosomal_contigs
      Array[String]? allosomal_contigs

      ############################
      #### Runtime attributes ####
      ############################
      RuntimeAttr? runtime_attr_ploidy
      RuntimeAttr? runtime_attr_case
      RuntimeAttr? runtime_attr_postprocess
      RuntimeAttr? runtime_attr_explode
    }

    call DetermineGermlineContigPloidyCaseMode {
        input:
            read_count_files = counts,
            contig_ploidy_model_tar = contig_ploidy_model_tar,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mapping_error_rate = ploidy_mapping_error_rate,
            sample_psi_scale = ploidy_sample_psi_scale,
            runtime_attr_override = runtime_attr_ploidy
    }

    scatter (scatter_index in range(length(gcnv_model_tars))) {
        call GermlineCNVCallerCaseMode {
            input:
                scatter_index = scatter_index,
                read_count_files = counts,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                gcnv_model_tar = gcnv_model_tars[scatter_index],
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                p_alt = gcnv_p_alt,
                cnv_coherence_length = gcnv_cnv_coherence_length,
                max_copy_number = gcnv_max_copy_number,
                mapping_error_rate = gcnv_mapping_error_rate,
                sample_psi_scale = gcnv_sample_psi_scale,
                depth_correction_tau = gcnv_depth_correction_tau,
                copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                learning_rate = gcnv_learning_rate,
                adamax_beta_1 = gcnv_adamax_beta_1,
                adamax_beta_2 = gcnv_adamax_beta_2,
                log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs = gcnv_min_training_epochs,
                max_training_epochs = gcnv_max_training_epochs,
                initial_temperature = gcnv_initial_temperature,
                num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                max_calling_iters = gcnv_max_calling_iters,
                caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                disable_annealing = gcnv_disable_annealing,
                runtime_attr_override = runtime_attr_case
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCaseMode.gcnv_call_tars)

    scatter (sample_index in range(length(counts))) {
        call CNVTasks.PostprocessGermlineCNVCalls {
            input:
                entity_id = count_entity_ids[sample_index],
                gcnv_calls_tars = call_tars_sample_by_shard[sample_index],
                gcnv_model_tars = gcnv_model_tars,
                calling_configs = GermlineCNVCallerCaseMode.calling_config_json,
                denoising_configs = GermlineCNVCallerCaseMode.denoising_config_json,
                gcnvkernel_version = GermlineCNVCallerCaseMode.gcnvkernel_version_json,
                sharded_interval_lists = GermlineCNVCallerCaseMode.sharded_interval_list,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                allosomal_contigs = allosomal_contigs,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                sample_index = sample_index,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                runtime_attr_override = runtime_attr_postprocess
        }
    }

    call CNVTasks.ExplodePloidyCalls {
        input :
            contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
            samples = count_entity_ids,
            linux_docker = linux_docker,
            runtime_attr_override = runtime_attr_explode
    }

    output {
        File contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar
        Array[File] sample_contig_ploidy_calls_tars = ExplodePloidyCalls.sample_contig_ploidy_calls_tar
        Array[Array[File]] gcnv_calls_tars = GermlineCNVCallerCaseMode.gcnv_call_tars
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCaseMode.gcnv_tracking_tar
        Array[File] genotyped_intervals_vcf = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        Array[File] genotyped_segments_vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf
        Array[File] denoised_copy_ratios = PostprocessGermlineCNVCalls.denoised_copy_ratios
    }
}

task DetermineGermlineContigPloidyCaseMode {
    input {
      Array[File] read_count_files
      File contig_ploidy_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale

      # Runtime parameters
      String gatk_docker
      RuntimeAttr? runtime_attr_override
    }

    Int base_disk_space_gb = 10
    Float unzip_factor = 10.0
    Int disk_gb = base_disk_space_gb + ceil(unzip_factor * (
        size(read_count_files, "GiB") + size(gatk4_jar_override, "GiB")
    ))

    RuntimeAttr default_attr = object {
      cpu_cores: 4,
      mem_gb: 8.5,
      disk_gb: disk_gb,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 * 0.6)
    Int cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -euo pipefail
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{cpu}
        export OMP_NUM_THREADS=~{cpu}

        mkdir input-contig-ploidy-model
        tar xzf ~{contig_ploidy_model_tar} -C input-contig-ploidy-model

        read_count_files_list=~{write_lines(read_count_files)}
        grep gz$ $read_count_files_list | xargs -l1 -P0 gunzip
        sed 's/\.gz$//' $read_count_files_list | \
            awk '{print "--input "$0}' > read_count_files.args

        gatk --java-options "-Xmx~{command_mem_mb}m" DetermineGermlineContigPloidy \
            --arguments_file read_count_files.args \
            --model input-contig-ploidy-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar c -C ~{output_dir_}/case-calls . | gzip -1 > case-contig-ploidy-calls.tar.gz
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File contig_ploidy_calls_tar = "case-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCaseMode {
    input {
      Int scatter_index
      Array[File] read_count_files
      File contig_ploidy_calls_tar
      File gcnv_model_tar
      String? output_dir
      File? gatk4_jar_override

      # Caller parameters
      Float? p_alt
      Float? cnv_coherence_length
      Int? max_copy_number

      # Denoising model parameters
      Float? mapping_error_rate
      Float? sample_psi_scale
      Float? depth_correction_tau
      String? copy_number_posterior_expectation_mode
      Int? active_class_padding_hybrid_mode

      # Hybrid ADVI parameters
      Float? learning_rate
      Float? adamax_beta_1
      Float? adamax_beta_2
      Int? log_emission_samples_per_round
      Float? log_emission_sampling_median_rel_error
      Int? log_emission_sampling_rounds
      Int? max_advi_iter_first_epoch
      Int? max_advi_iter_subsequent_epochs
      Int? min_training_epochs
      Int? max_training_epochs
      Float? initial_temperature
      Int? num_thermal_advi_iters
      Int? convergence_snr_averaging_window
      Float? convergence_snr_trigger_threshold
      Int? convergence_snr_countdown_window
      Int? max_calling_iters
      Float? caller_update_convergence_threshold
      Float? caller_internal_admixing_rate
      Float? caller_external_admixing_rate
      Boolean? disable_annealing

      # Runtime parameters
      String gatk_docker
      RuntimeAttr? runtime_attr_override
    }

    Int base_disk_space_gb = 10
    Float unzip_factor = 10.0
    Int disk_gb = base_disk_space_gb + ceil(unzip_factor * (
        size(read_count_files, "GiB") + size(contig_ploidy_calls_tar, "GiB") +
        size(gcnv_model_tar, "GiB") + size(gatk4_jar_override, "GiB")
    ))

    RuntimeAttr default_attr = object {
      cpu_cores: 4,
      mem_gb: 10,
      disk_gb: disk_gb,
      boot_disk_gb: 10,
      preemptible_tries: 3,
      max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int command_mem_mb = ceil(mem_gb * 1000 * 0.6)
    Int cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -euo pipefail
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{cpu}
        export OMP_NUM_THREADS=~{cpu}

        mkdir contig-ploidy-calls-dir
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls-dir

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        read_count_files_list=~{write_lines(read_count_files)}
        grep gz$ "$read_count_files_list" | xargs -l1 -P0 gunzip
        sed 's/\.gz$//' "$read_count_files_list" \
            | awk '{print "--input "$0}' \
            > read_count_files.args

        function get_seeded_random() {
            SEED="$1"
            openssl enc -aes-256-ctr -pass pass:"$SEED" -nosalt </dev/zero 2>/dev/null
        }

        function run_gcnv_case() {
            gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
                --run-mode CASE \
                --arguments_file read_count_files.args \
                --contig-ploidy-calls contig-ploidy-calls-dir \
                --model gcnv-model \
                --output ~{output_dir_} \
                --output-prefix case \
                --verbosity DEBUG \
                --p-alt ~{default="1e-6" p_alt} \
                --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
                --max-copy-number ~{default="5" max_copy_number} \
                --mapping-error-rate ~{default="0.01" mapping_error_rate} \
                --sample-psi-scale ~{default="0.0001" sample_psi_scale} \
                --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
                --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
                --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
                --learning-rate ~{default="0.05" learning_rate} \
                --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
                --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
                --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
                --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
                --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
                --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
                --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
                --min-training-epochs ~{default="10" min_training_epochs} \
                --max-training-epochs ~{default="100" max_training_epochs} \
                --initial-temperature ~{default="2.0" initial_temperature} \
                --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
                --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
                --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
                --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
                --max-calling-iters ~{default="10" max_calling_iters} \
                --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
                --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
                --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
                --disable-annealing ~{default="false" disable_annealing}
        }

        {
            # Try to run gcnv case mode. Rarely, a bad initial seed results in NaN errors...
            run_gcnv_case
        } || {
            # shuffle input arguments in a deterministic manner, resulting in a new seed
            shuf --random-source=<(get_seeded_random 42) --output=read_count_files.args read_count_files.args
            # run gcnv case mode one more time
            run_gcnv_case
        }

        tar c -C ~{output_dir_}/case-tracking . | gzip -1 > case-gcnv-tracking-~{scatter_index}.tar.gz

        # tar output calls, ensuring output files are numbered in the same order as original sample list
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        CURRENT_SAMPLE=0
        sed 's/\.gz$//' "$read_count_files_list" \
            | while read READ_COUNT_FILE; do
                SAMPLE_NAME=$(zgrep "^@RG" "$READ_COUNT_FILE" | cut -d: -f3)
                SAMPLE_PATH=$(dirname $(grep -lR -m1 "^$SAMPLE_NAME$" "~{output_dir_}/case-calls"))
                CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
                tar czf case-gcnv-calls-shard-~{scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz \
                    -C "$SAMPLE_PATH" .
                ((++CURRENT_SAMPLE))
            done
    >>>
    runtime {
      cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
      memory: mem_gb + " GiB"
      disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
      docker: gatk_docker
      preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        Array[File] gcnv_call_tars = glob("case-gcnv-calls-shard-~{scatter_index}-sample-*.tar.gz")
        File gcnv_tracking_tar = "case-gcnv-tracking-~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"
    }
}
