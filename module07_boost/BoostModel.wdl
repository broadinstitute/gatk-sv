##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks
import "SplitPerSiteVCF.wdl" as split_per_site_vcf
import "SplitPerSampleGTGQandBEDPerSampleList.wdl" as split_per_sample_gtgq_and_bed
import "AnnotateILFeaturesBed.wdl" as annotate_il_features
import "AnnotateTrainingFeatures.wdl" as annotate_training_features

workflow BoostModel{
    input{
        Array[File] cleanVcfs
        Array[File] cleanVcfIdxes
        Array[File]? cleanBeds
        Array[String] VcfPrefixes

        Array[File] cluster_manta_files
        Array[File] cluster_wham_files
        Array[File] cluster_melt_files
        Array[File] cluster_depth_files
        Array[String] cluster_names

        Array[String] training_samples
        Array[File] vapor_result_files

        File simp_rep
        File seg_dup
        File rep_mask

        String prefix

        File ref_fasta
        File ref_fai
        File ref_dict
        Int min_shard_size
        Array[String] contigs

        String sv_types
        String size_ranges
        String af_ranges

        String vapor_docker
        String duphold_docker
        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        Boolean requester_pays_crams = false
        Boolean run_genomic_context_anno = false
        Boolean run_extract_algo_evi = false
        Boolean run_duphold = false
        Boolean run_extract_gt_gq = true
        Boolean run_versus_raw_vcf = true
        Boolean run_rdpesr_anno = false
        Boolean run_vapor = false

        RuntimeAttr? runtime_attr_vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_concat_refs
        RuntimeAttr? runtime_attr_concat_gtgqs
        RuntimeAttr? runtime_attr_override_vcf2bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
        RuntimeAttr? runtime_split_per_sample_gtgq
        RuntimeAttr? runtime_attr_override_format_ref 
        RuntimeAttr? runtime_attr_override_concat_beds
        RuntimeAttr? runtime_attr_override_inte_features
        RuntimeAttr? runtime_attr_override_extract_sample_list
        RuntimeAttr? runtime_attr_bed_vs_ont
        RuntimeAttr? runtime_attr_bed_vs_hgsv
        RuntimeAttr? runtime_attr_bed_vs_array
        RuntimeAttr? runtime_attr_bed_vs_pacbio
        RuntimeAttr? runtime_attr_bed_vs_bionano
        RuntimeAttr? runtime_attr_override_inte_anno
        RuntimeAttr? runtime_attr_override_extract_training_samples
        RuntimeAttr? runtime_attr_override_organize_training_data
        RuntimeAttr? runtime_attr_override_organize_training_data_per_sample
        RuntimeAttr? runtime_attr_override_train_boost_model
        RuntimeAttr? runtime_attr_override_apply_boost_model
    }

    call split_per_site_vcf.SplitPerSiteVCF as SplitPerSiteVCF{
        input:
            cleanVcfs = cleanVcfs, 
            simp_rep = simp_rep,
            seg_dup = seg_dup,
            rep_mask = rep_mask,
            prefix = prefix,

            sv_pipeline_docker = sv_pipeline_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_benchmark_docker = sv_benchmark_docker,

            runtime_attr_override_vcf2bed = runtime_attr_override_vcf2bed,
            runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
            runtime_attr_override_inte_gc = runtime_attr_override_inte_gc,
            runtime_attr_override_concat_beds = runtime_attr_override_concat_beds,
            runtime_attr_override_inte_features = runtime_attr_override_inte_features
    }

    scatter(i in range(length(cluster_manta_files))){
        call ExtractSampleAndVcfList{
            input:
                cluster_manta = cluster_manta_files[i],
                cluster_wham = cluster_wham_files[i],
                cluster_melt = cluster_melt_files[i],
                cluster_depth = cluster_depth_files[i],
                prefix = prefix,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_extract_sample_list
        }

        call split_per_sample_gtgq_and_bed.SplitPerSampleGTGQandBEDPerSampleList as SplitPerSampleGTGQandBEDPerSampleList{
            input:
                cleanVcfs = cleanVcfs,
                cleanVcfIdxes = cleanVcfIdxes,
                cleanBeds = cleanBeds,
                SampleList = ExtractSampleAndVcfList.sample_list,
                prefixes = VcfPrefixes,

                sv_benchmark_docker = sv_benchmark_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_split_per_sample_gtgq = runtime_split_per_sample_gtgq,
                runtime_attr_concat_gtgqs = runtime_attr_concat_gtgqs,
                runtime_attr_concat_refs = runtime_attr_concat_refs,
                runtime_attr_override_vcf2bed = runtime_attr_override_vcf2bed
        }

        Array[String] sample_list_extracted = read_lines(ExtractSampleAndVcfList.sample_list)

        call annotate_il_features.AnnotateILFeatures as AnnotateILFeatures{
            input:
                samples = sample_list_extracted,
                raw_mantas = ExtractSampleAndVcfList.manta_list,
                raw_whams = ExtractSampleAndVcfList.wham_list,
                raw_melts = ExtractSampleAndVcfList.melt_list,
                raw_depths = ExtractSampleAndVcfList.depth_list,
                gtgqs = SplitPerSampleGTGQandBEDPerSampleList.gtgq_out,
                beds = SplitPerSampleGTGQandBEDPerSampleList.bed_out,
                prefix = cluster_names[i],
                sv_benchmark_docker = sv_benchmark_docker,
                duphold_docker = duphold_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    call ExtractTrainingSamples{
        input:
            training_samples = training_samples,
            il_anno_tarballs = AnnotateILFeatures.anno_tar,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_extract_training_samples
    }

    call annotate_training_features.AnnotateTrainingFeatures as AnnotateTrainingFeatures{
        input:
            samples = training_samples,
            IL_anno_files = ExtractTrainingSamples.training_anno,
            
            vapor_result_files = vapor_result_files, 

            run_vapor = run_vapor,
            run_duphold = run_duphold,
            run_rdpesr_anno = run_rdpesr_anno ,
            run_extract_gt_gq = run_extract_gt_gq,
            run_versus_raw_vcf = run_versus_raw_vcf,
            requester_pays_crams = requester_pays_crams,
            run_extract_algo_evi = run_extract_algo_evi,
            run_genomic_context_anno = run_genomic_context_anno,

            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            min_shard_size = min_shard_size,

            vapor_docker = vapor_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            sv_benchmark_docker = sv_benchmark_docker,

            runtime_attr_vapor = runtime_attr_vapor,
            runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
            runtime_attr_vcf2bed = runtime_attr_vcf2bed,
            runtime_attr_SplitVcf = runtime_attr_SplitVcf,
            runtime_attr_ConcatBeds = runtime_attr_ConcatBeds,
            runtime_attr_bed_vs_ont = runtime_attr_bed_vs_ont,
            runtime_attr_bed_vs_hgsv = runtime_attr_bed_vs_hgsv,
            runtime_attr_bed_vs_array = runtime_attr_bed_vs_array,
            runtime_attr_bed_vs_pacbio = runtime_attr_bed_vs_pacbio,
            runtime_attr_bed_vs_bionano = runtime_attr_bed_vs_bionano,
            runtime_attr_override_inte_anno = runtime_attr_override_inte_anno,
            runtime_attr_override_format_ref = runtime_attr_override_format_ref
    }

    scatter(i in range(length(training_samples))){
        call OrganizeTrainingDataPerSample{
            input:
                sample = training_samples[i],
                site_anno = SplitPerSiteVCF.anno_bed,
                training_anno = AnnotateTrainingFeatures.annotated_files[i],
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_organize_training_data_per_sample
        }
    }

    call OrganizeTrainingData{
        input:
            training_per_sample = OrganizeTrainingDataPerSample.organized_training_anno,
            prefix = prefix,
            sv_types = sv_types,
            size_ranges = size_ranges,
            af_ranges = af_ranges,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_organize_training_data
    }

    scatter(train_data in OrganizeTrainingData.training_data){
        call TrainBoostModel{
            input:
                training_data = train_data,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_train_boost_model

        }
    }

    scatter(il_anno_list in AnnotateILFeatures.anno_tar){
        call TestBoostModel{
            input:
                trained_models = TrainBoostModel.trained_model,
                test_data = il_anno_list,
                site_anno = SplitPerSiteVCF.anno_bed,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_apply_boost_model
        }
    }

    output{
        Array[File] bs_models = TrainBoostModel.trained_model
        Array[File] bs_filtered_tarballs = TestBoostModel.bs_filtered_tar
    }

}



task SplitSampleList{
    input{
        File manifest
        String prefix
        Int count
        String sv_base_mini_docker
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


    command<<<
        split -l ~{count} -a 6 <(cut -f1 ~{manifest}) ~{prefix}.
        split -l ~{count} -a 6 <(cut -f2 ~{manifest}) ~{prefix}.manta.
        split -l ~{count} -a 6 <(cut -f3 ~{manifest}) ~{prefix}.wham.
        split -l ~{count} -a 6 <(cut -f4 ~{manifest}) ~{prefix}.melt.

    >>>

    output{
        Array[File] sample_list = glob("~{prefix}.*")
        Array[File] manta_list = glob("~{prefix}.manta.*")
        Array[File] wham_list = glob("~{prefix}.wham.*")
        Array[File] melt_list = glob("~{prefix}.melt.*")
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractSampleAndVcfList{
    input{
        File cluster_manta
        File cluster_wham
        File cluster_melt
        File cluster_depth
        String prefix
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


    command<<<
        set -eux

        zcat ~{cluster_manta} | grep -m1 CHROM | cut -f10- | sed -e 's/\t/\n/g' > sample_list.tsv

        set -o pipefail

        svtk vcf2bed -i SVTYPE -i SVLEN ~{cluster_manta} ~{prefix}.manta.bed
        svtk vcf2bed -i SVTYPE -i SVLEN ~{cluster_wham} ~{prefix}.wham.bed
        svtk vcf2bed -i SVTYPE -i SVLEN ~{cluster_melt} ~{prefix}.melt.bed
        svtk vcf2bed -i SVTYPE -i SVLEN ~{cluster_depth} ~{prefix}.depth.bed
        head -1 ~{prefix}.manta.bed | cut -f1-4,7,8 > header_query

        i=0
        while read sample;
        do
            app=$(printf "%06d" $i)
            cat header_query <(grep $sample ~{prefix}.manta.bed | cut -f1-4,7,8) | bgzip > ~{prefix}.manta.$app.$sample.query.gz
            cat header_query <(grep $sample ~{prefix}.wham.bed  | cut -f1-4,7,8) | bgzip > ~{prefix}.wham.$app.$sample.query.gz
            cat header_query <(grep $sample ~{prefix}.melt.bed  | cut -f1-4,7,8) | bgzip > ~{prefix}.melt.$app.$sample.query.gz
            cat header_query <(grep $sample ~{prefix}.depth.bed | cut -f1-4,7,8) | bgzip > ~{prefix}.depth.$app.$sample.query.gz
            i=$((i + 1))
        done < sample_list.tsv
    >>>

    output{
        File sample_list = "sample_list.tsv"
        Array[File] manta_list = glob("~{prefix}.manta.*")
        Array[File] wham_list = glob("~{prefix}.wham.*")
        Array[File] melt_list = glob("~{prefix}.melt.*")
        Array[File] depth_list = glob("~{prefix}.depth.*")
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

task ExtractTrainingSamples{
    input{
        Array[String] training_samples
        Array[File] il_anno_tarballs
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


    command<<<
        set -Eeuo pipefail
        mkdir tmp/
        while read SPLIT; do
            tar zxvf $SPLIT
            mv merged/* tmp/
        done < ~{write_lines(il_anno_tarballs)}

        mkdir target/
        i=0
        while read SAMPLE; do
            app=$(printf "%06d" $i)
            mv tmp/*$SAMPLE*  target/$app.$SAMPLE.anno.bed.gz
            i=$((i + 1))
        done < ~{write_lines(training_samples)}

    >>>

    output{
        Array[File] training_anno = glob("target/*.anno.bed.gz")
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

task OrganizeTrainingDataPerSample{
    input{
        String sample
        File site_anno
        File training_anno
        String sv_benchmark_docker
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


    command<<<
        set -Eeuo pipefail

        Rscript /src/integrate_ref_support.R \
        --site_anno ~{site_anno} \
        --sample_anno ~{training_anno} \
        --output  ~{sample}.training
        bgzip ~{sample}.training
    >>>

    output{
        File organized_training_anno = "~{sample}.training.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task OrganizeTrainingData{
    input{
        Array[File] training_per_sample
        String prefix
        String sv_types
        String size_ranges
        String af_ranges
        String sv_benchmark_docker
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


    command<<<
        set -eux

        # note head -n1 stops reading early and sends SIGPIPE to zcat,
        # so setting pipefail here would result in early termination
        
        zcat ~{training_per_sample[0]} | head -n1 > integrated_training.tsv

        # no more early stopping
        set -o pipefail

        while read SPLIT; do
          zcat $SPLIT | tail -n+2
        done < ~{write_lines(training_per_sample)} \
          >> integrated_training.tsv

        python /src/split_ref_support.py \
            integrated_training.tsv \
            ~{prefix} \
            -t ~{sv_types} \
            -s ~{size_ranges} \
            -a ~{af_ranges}

        >>>

    output{
        File integrated_training = "integrated_training.tsv"
        Array[File] training_data = glob("~{prefix}.*.training.gz")
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task TrainBoostModel{
    input{
        File training_data
        String sv_benchmark_docker
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

    String prefix = basename(training_data, '.training.gz')
    command<<<
        set -Eeuo pipefail
        Rscript /src/TrainRandomForestModel.R \
            --train_data ~{training_data} \
            --output ~{prefix}.RDS

        >>>

    output{
        File trained_model = "~{prefix}.RDS"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task TestBoostModel{
    input{
        Array[File] trained_models
        File test_data
        File site_anno
        String sv_benchmark_docker
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

    String prefix = basename(test_data, '.tar.gz')
    command<<<
        set -Eeuo pipefail

        tar zxvf ~{test_data}
        ls merged/*anno.bed.gz > test_data_list.tsv

        Rscript /src/apply_boost_models.R \
            --input test_data_list.tsv \
            --site ~{site_anno} \
            --Boost_model_del_sz1_af1 ~{trained_models[0]} \
            --Boost_model_del_sz1_af2 ~{trained_models[1]} \
            --Boost_model_del_sz1_af3 ~{trained_models[2]} \
            --Boost_model_del_sz1_af4 ~{trained_models[3]} \
            --Boost_model_del_sz2_af1 ~{trained_models[4]} \
            --Boost_model_del_sz2_af2 ~{trained_models[5]} \
            --Boost_model_del_sz2_af3 ~{trained_models[6]} \
            --Boost_model_del_sz2_af4 ~{trained_models[7]} \
            --Boost_model_del_sz3_af1 ~{trained_models[8]} \
            --Boost_model_del_sz3_af2 ~{trained_models[9]} \
            --Boost_model_del_sz3_af3 ~{trained_models[10]} \
            --Boost_model_del_sz3_af4 ~{trained_models[11]} \
            --Boost_model_del_sz4_af1 ~{trained_models[12]} \
            --Boost_model_del_sz4_af2 ~{trained_models[13]} \
            --Boost_model_del_sz4_af3 ~{trained_models[14]} \
            --Boost_model_del_sz4_af4 ~{trained_models[15]} \
            --Boost_model_dup_sz1_af1 ~{trained_models[16]} \
            --Boost_model_dup_sz1_af2 ~{trained_models[17]} \
            --Boost_model_dup_sz1_af3 ~{trained_models[18]} \
            --Boost_model_dup_sz1_af4 ~{trained_models[19]} \
            --Boost_model_dup_sz2_af1 ~{trained_models[20]} \
            --Boost_model_dup_sz2_af2 ~{trained_models[21]} \
            --Boost_model_dup_sz2_af3 ~{trained_models[22]} \
            --Boost_model_dup_sz2_af4 ~{trained_models[23]} \
            --Boost_model_dup_sz3_af1 ~{trained_models[24]} \
            --Boost_model_dup_sz3_af2 ~{trained_models[25]} \
            --Boost_model_dup_sz3_af3 ~{trained_models[26]} \
            --Boost_model_dup_sz3_af4 ~{trained_models[27]} \
            --Boost_model_dup_sz4_af1 ~{trained_models[28]} \
            --Boost_model_dup_sz4_af2 ~{trained_models[29]} \
            --Boost_model_dup_sz4_af3 ~{trained_models[30]} \
            --Boost_model_dup_sz4_af4 ~{trained_models[31]} \
            --Boost_model_ins_sz1 ~{trained_models[32]} \
            --Boost_model_ins_sz2 ~{trained_models[33]} \
            --Boost_model_ins_sz3 ~{trained_models[34]} \
            --Boost_model_ins_sz4 ~{trained_models[35]} \
            --Boost_model_ins_alu ~{trained_models[36]} \
            --Boost_model_ins_l1 ~{trained_models[37]} \
            --Boost_model_ins_sva ~{trained_models[38]} \
            --Boost_model_ins_me ~{trained_models[39]} \
            --Boost_model_inv ~{trained_models[40]}

        mkdir bs_filtered/
        mv merged/*.bs bs_filtered/
        tar -czvf ~{prefix}.bs_filtered.tar.gz bs_filtered/
        >>>

    output{
        File bs_filtered_tar = "~{prefix}.bs_filtered.tar.gz"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
