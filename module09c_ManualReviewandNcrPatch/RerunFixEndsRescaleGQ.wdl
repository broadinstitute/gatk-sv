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

workflow RerunFixEndsRescaleGQ{
    input{ 
        File vcf 
        File vcf_idx

        String prefix
        String i

        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_fix_bad_ends 
        }
 
    call MiniTasks.FixEndsRescaleGQ{
        input:
            vcf = vcf,
            prefix = "~{prefix}.~{i}",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_fix_bad_ends
    }
}




