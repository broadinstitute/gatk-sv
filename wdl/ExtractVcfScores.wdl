version 1.0

import "Utils.wdl" as Utils
import "TasksMakeCohortVcf.wdl" as TasksMakeCohortVcf
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties

workflow ExtractVcfScores {
    input {
        File vcf
        File vcf_index
        Array[String] wanted_properties
        Boolean collect_whole_vcf_parquet = true
        Boolean compute_weights = false
        Array[String] tarred_properties_to_parquet_args = []
        String samtools_cloud_docker
        String sv_utils_docker
        Int records_per_shard = 10000
    }

    call TasksMakeCohortVcf.ScatterVcf {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            prefix=sub(sub(basename(vcf), ".gz$", ""), ".vcf$", ""),
            records_per_shard=records_per_shard,
            sv_pipeline_docker=samtools_cloud_docker
    }

    call Utils.GetVcfSize {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            samtools_cloud_docker=samtools_cloud_docker
    }

    scatter(vcf_shard in ScatterVcf.shards) {
        call Utils.IndexVcf {
            input:
                vcf=vcf_shard,
                samtools_cloud_docker=samtools_cloud_docker
        }
        call VcfToParquet {
            input:
                vcf=vcf_shard,
                wanted_properties=wanted_properties,
                sv_utils_docker=sv_utils_docker
        }

    }

    if(collect_whole_vcf_parquet) {
        # get parquet file for whole data set
        call ExtractVcfGqFilterProperties.ConcatenateParquetTars {
            input:
                parquet_tars=VcfToParquet.properties_parquet_tar,
                output_base_name=sub(sub(basename(vcf), ".gz$", ""), ".vcf$", ""),
                compute_weights=compute_weights,
                sv_utils_docker=sv_utils_docker
        }
    }

    output {
        Array[File] vcf_shards = ScatterVcf.shards
        Array[File] properties_parquet_shards = VcfToParquet.properties_parquet_tar
        File? properties_parquet_tar = ConcatenateParquetTars.combined_parquet_tar
        Int num_samples = GetVcfSize.num_samples
    }
}

task VcfToParquet {
    input {
        File vcf
        Array[String] wanted_properties
        String sv_utils_docker
        Float mem_gb = 8
    }

    Int disk_gb = round(50 + size([vcf], "GiB"))

    # create output filename:
    #   strip .vcf.gz from end (if present)
    #   add wanted properties (separated by underscores)
    #   and add ".pickle.bz2"
    #   (sep doesn't work the way it's supposed to outside of command block, so move part of this logic to command)
    String vcf_basename = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "")
    String parquet_file_name = vcf_basename + ".pq.tar"

    runtime {
        docker: sv_utils_docker
        cpu: 8
        preemptible: 1
        max_retries: 1
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
    gq-recalibrator vcf-to-dask-parquet \
        -i ~{vcf} \
        -o ~{parquet_file_name} \
        --wanted-properties "~{sep="," wanted_properties}" \
        --drop-rows-with-missing-values \
        --temp-dir "$(mktemp -d --tmpdir=.)"

    >>>

    output {
        File properties_parquet_tar = parquet_file_name
    }
}