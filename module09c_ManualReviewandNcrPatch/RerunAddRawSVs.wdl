##########################################################################################

## Copyright Broad Institute, 2022
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
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results
import "ShardedManualReview.wdl" as sharded_manual_review

workflow RerunAddRawSVs{
    input{ 
        File reviewed_vcf 
        File reviewed_vcf_idx
        File raw_SVs

        String prefix
        String contig

        String sv_base_mini_docker
        String sv_benchmark_docker

        RuntimeAttr? runtime_attr_add_raw_SVs 
        RuntimeAttr? runtime_attr_split_raw_SVs_per_chr
        }
 

    if (contig!="whole_genome"){
        call sharded_manual_review.SplitRawSVsPerChr{
            input:
                raw_SVs = raw_SVs,
                contig = contig,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_raw_SVs_per_chr
        }
    }

    File split_raw_SVs = select_first([SplitRawSVsPerChr.raw_SV_per_chr, raw_SVs])

    call revise_vcf_with_manual_results.AddRawSVs{
        input:
            prefix = prefix,
            batch_name = contig,
            vcf_file = reviewed_vcf,
            raw_SVs = split_raw_SVs,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_add_raw_SVs
    }

    output{
        File revised_output_vcf = AddRawSVs.vcf_with_raw_SVs
        File revised_output_vcf_idx = AddRawSVs.vcf_idx_with_raw_SVs

    }
}

