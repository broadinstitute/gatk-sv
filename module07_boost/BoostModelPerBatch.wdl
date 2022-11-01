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

workflow BoostModelPerBatch{
    input{
        Array[File] cleanVcfs
        Array[File] cleanVcfIdxes
        Array[File]? cleanBeds
        Array[String] VcfPrefixes

        File cluster_manta
        File cluster_wham
        File cluster_melt
        File cluster_depth
        String cluster_prefix

        File simp_rep
        File seg_dup
        File rep_mask

        Int size_sample_list
        String prefix

        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String duphold_docker

        RuntimeAttr? runtime_attr_override_vcf2bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
        RuntimeAttr? runtime_split_per_sample_gtgq
        RuntimeAttr? runtime_attr_concat_gtgqs
        RuntimeAttr? runtime_attr_concat_refs
        RuntimeAttr? runtime_attr_override_inte_features
        RuntimeAttr? runtime_attr_override_concat_beds
        RuntimeAttr? runtime_attr_override_extract_sample_list

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

    call ExtractSampleAndVcfList{
        input:
            cluster_manta = cluster_manta,
            cluster_wham = cluster_wham,
            cluster_melt = cluster_melt,
            cluster_depth = cluster_depth,
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
            prefix = cluster_prefix,
            samples = sample_list_extracted,
            raw_mantas = ExtractSampleAndVcfList.manta_list,
            raw_whams = ExtractSampleAndVcfList.wham_list,
            raw_melts = ExtractSampleAndVcfList.melt_list,
            raw_depths = ExtractSampleAndVcfList.depth_list,
            gtgqs = SplitPerSampleGTGQandBEDPerSampleList.gtgq_out,
            beds = SplitPerSampleGTGQandBEDPerSampleList.bed_out,
            sv_benchmark_docker = sv_benchmark_docker,
            duphold_docker = duphold_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker
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


