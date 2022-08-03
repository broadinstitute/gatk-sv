version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"

# Workflow to annotate vcf file with genomic context
workflow ApplyRandomForestModels {
    input {
        Array[File] PBSV_RF_models # RF_models: [DEL_under100bp,DEL_100to1Kb,DEL_1to5Kb,DEL_over5Kb,DUP_under100bp,DUP_100to1Kb,DUP_1to5Kb,DUP_over5Kb,INS_under100bp,INS_100to1Kb,INS_1to5Kb,INS_over5Kb,ALU,L1,SVA]
        Array[File] VaPoR_RF_models

        Array[File] site_level_annotations #[gnomad-sv-v3.final_cleanup.anno.rds, gnomad-sv-v3.final_cleanup.info4.rds, gnomAD_SV_v3.posthoc_filter_integration.nonredundant_union.integrated.pb_samples.SVID_NCR.rds ]

        Array[File] remove_redun_annotations #[ ../../INS_recluster/Redun_Removed.SR.INS.SVID_anno.V2, ../../INS_recluster/Redun_Removed.SR.INS.SVID_vs_closest_MEMBER.tsv, ../../INS_recluster/Redun_Removed.SR.INS.sites]

        Array[String] sample_list
        Array[File] loose_union
        Array[File] GQ_reacali
        Array[File] boost_dyna
        Array[File] boost_fix
        Array[File] minGQ10
        Array[File] boost_score

        String sv_base_mini_docker
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override_integrate_test_data
        RuntimeAttr? runtime_attr_override_apply_random_forest_model
        }

    call IntegrateTestingData{
        input:
            prefix = "Integration",
            sample_list = sample_list,
            loose_union = loose_union,
            GQ_reacali = GQ_reacali,
            boost_dyna = boost_dyna,
            boost_fix = boost_fix,
            minGQ10 = minGQ10,
            boost_score = boost_score,
            site_level_annotations = site_level_annotations,
            remove_redun_annotations = remove_redun_annotations,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_integrate_test_data    
    }

    call ApplyRandomForestModels{
        input:
            prefix = "RF_results",
            PBSV_RF_models = PBSV_RF_models, 
            VaPoR_RF_models = VaPoR_RF_models,
            Input_data = IntegrateTestingData.test_dat,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_apply_random_forest_model
    }

    output{
        Array[File] integrated_test = IntegrateTestingData.test_dat
        Array[File] random_forest_output = ApplyRandomForestModels.RF_results
    }
}

task IntegrateTestingData {
    input {
        Array[String] sample_list
        Array[File] loose_union
        Array[File] GQ_reacali
        Array[File] boost_dyna
        Array[File] boost_fix
        Array[File] minGQ10
        Array[File] boost_score
        Array[File] site_level_annotations
        Array[File] remove_redun_annotations
        String prefix
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 7.5,
        disk_gb: 10,
        cpu_cores: 1,
        preemptible_tries: 1,
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
        docker: sv_benchmark_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        while read name; do
            echo $name
        done < ~{write_lines(sample_list)} > sample_list

        while read name; do
            echo $name
        done < ~{write_lines(loose_union)} > loose_union_list

        while read name; do
            echo $name
        done < ~{write_lines(GQ_reacali)} > GQ_reacali_list

        while read name; do
            echo $name
        done < ~{write_lines(boost_dyna)} > boost_dyna_list

        while read name; do
            echo $name
        done < ~{write_lines(boost_fix)} > boost_fix_list

        while read name; do
            echo $name
        done < ~{write_lines(minGQ10)} > minGQ10_list

        while read name; do
            echo $name
        done < ~{write_lines(boost_score)} > boost_score_list

        while read name; do
            echo ~{prefix}.${name}.test_data.tsv
        done < ~{write_lines(sample_list)} > output_list

        Rscript /src/integrate_per_sample_testdata.test.R \
            --sample sample_list \
            --loose_union loose_union_list \
            --GQ_reacali GQ_reacali_list \
            --boost_dyna boost_dyna_list \
            --boost_fix boost_fix_list \
            --minGQ10perc minGQ10_list \
            --boost_score boost_score_list \
            --output  output_list \
            --SVID_anno ~{site_level_annotations[0]} \
            --SVID_coor ~{site_level_annotations[1]} \
            --SVID_NCR  ~{site_level_annotations[2]} \
            --ReCluster_anno ~{remove_redun_annotations[0]} \
            --ReCluster_closest ~{remove_redun_annotations[1]} \
            --ReCluster_sites ~{remove_redun_annotations[2]}

    >>>

    output {
        Array[File] test_dat = glob("~{prefix}.*")
    }
}

task ApplyRandomForestModels{
    input{
        Array[File] PBSV_RF_models
        Array[File] VaPoR_RF_models
        Array[File] Input_data
        String prefix
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 5,
        disk_gb: 10,
        cpu_cores: 1,
        preemptible_tries: 1,
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
        docker: sv_benchmark_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        while read name; do
            echo $name
        done < ~{write_lines(Input_data)} > input_data_list

        Rscript /src/apply_randon_forest_models.R \
            --input input_data_list \
            --PBSV_model_del_sz1   ~{PBSV_RF_models[0]} \
            --PBSV_model_del_sz2   ~{PBSV_RF_models[1]} \
            --PBSV_model_del_sz3   ~{PBSV_RF_models[2]} \
            --PBSV_model_del_sz4   ~{PBSV_RF_models[3]} \
            --PBSV_model_dup_sz1   ~{PBSV_RF_models[4]} \
            --PBSV_model_dup_sz2   ~{PBSV_RF_models[5]} \
            --PBSV_model_dup_sz3   ~{PBSV_RF_models[6]} \
            --PBSV_model_dup_sz4   ~{PBSV_RF_models[7]} \
            --PBSV_model_ins_sz1   ~{PBSV_RF_models[8]} \
            --PBSV_model_ins_sz2   ~{PBSV_RF_models[9]} \
            --PBSV_model_ins_sz3   ~{PBSV_RF_models[10]} \
            --PBSV_model_ins_sz4   ~{PBSV_RF_models[11]} \
            --PBSV_model_ins_alu   ~{PBSV_RF_models[12]} \
            --PBSV_model_ins_l1    ~{PBSV_RF_models[13]} \
            --PBSV_model_ins_sva   ~{PBSV_RF_models[14]} \
            --VaPoR_model_del_sz1  ~{VaPoR_RF_models[0]} \
            --VaPoR_model_del_sz2  ~{VaPoR_RF_models[1]} \
            --VaPoR_model_del_sz3  ~{VaPoR_RF_models[2]} \
            --VaPoR_model_del_sz4  ~{VaPoR_RF_models[3]} \
            --VaPoR_model_dup_sz1  ~{VaPoR_RF_models[4]} \
            --VaPoR_model_dup_sz2  ~{VaPoR_RF_models[5]} \
            --VaPoR_model_dup_sz3  ~{VaPoR_RF_models[6]} \
            --VaPoR_model_dup_sz4  ~{VaPoR_RF_models[7]} \
            --VaPoR_model_ins_sz1  ~{VaPoR_RF_models[8]} \
            --VaPoR_model_ins_sz2  ~{VaPoR_RF_models[9]} \
            --VaPoR_model_ins_sz3  ~{VaPoR_RF_models[10]} \
            --VaPoR_model_ins_sz4  ~{VaPoR_RF_models[11]} \
            --VaPoR_model_ins_alu  ~{VaPoR_RF_models[12]} \
            --VaPoR_model_ins_l1   ~{VaPoR_RF_models[13]} \
            --VaPoR_model_ins_sva  ~{VaPoR_RF_models[14]}
    >>>

    output{
        Array[File] RF_results = glob("~{prefix}.*")
    }

}



 
