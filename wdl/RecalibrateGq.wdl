version 1.0

import "RecalibrateGqShards.wdl" as RecalibrateGqShards
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties

workflow RecalibrateGq {
    input {
        File vcf
        File vcf_index
        File gq_recalibrator_model
        # recalibrated_vcf_basename must be a string with no .vcf.gz suffix
        String recalibrated_vcf_basename = sub(
            sub(basename(vcf), ".gz$", ""),
            ".vcf$", "-recalibrated"
        )
        Array[File] genome_tracks
        Array[String] recalibrate_gq_args = []
        Array[String] annotate_gq_args = []
        String sv_base_mini_docker
        String samtools_cloud_docker
        String sv_utils_docker
        String gatk_recalibrator_docker
    }


    call ExtractVcfGqFilterProperties.ExtractVcfGqFilterProperties {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            genome_tracks=genome_tracks,
            samtools_cloud_docker=samtools_cloud_docker,
            gatk_recalibrator_docker=gatk_recalibrator_docker,
            sv_utils_docker=sv_utils_docker
    }

    call RecalibrateGqShards.RecalibrateGqShards {
        input:
            vcf_shards=ExtractVcfGqFilterProperties.vcf_shards,
            properties_parquet_shards=ExtractVcfGqFilterProperties.properties_parquet_shards,
            gq_recalibrator_model=gq_recalibrator_model,
            recalibrated_vcf_basename=recalibrated_vcf_basename,
            num_samples=ExtractVcfGqFilterProperties.num_samples,
            recalibrate_gq_args=recalibrate_gq_args,
            annotate_gq_args=annotate_gq_args,
            sv_base_mini_docker=sv_base_mini_docker,
            sv_utils_docker=sv_utils_docker
    }

    output {
        File recalibrated_scores_parquet = RecalibrateGqShards.recalibrated_scores_parquet
        File recalibrated_vcf = RecalibrateGqShards.recalibrated_vcf
        File recalibrated_vcf_index = RecalibrateGqShards.recalibrated_vcf_index
    }
}
