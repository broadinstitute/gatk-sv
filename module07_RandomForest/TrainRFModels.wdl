version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "TrainRFModelsBySVtype.wdl" as TrainRFModelsBySVtype

# Workflow to annotate vcf file with genomic context
workflow TrainRandomForestModels {
    input {
        File training_data
        File inheritance
        File genomic_context

        Array[String] svtype_list
        Array[String] size_range_list

        String prefix

        String sv_benchmark_docker

        RuntimeAttr? runtime_attr_override_train_RF_model

    }

    scatter (svtype in svtype_list){
            call TrainRFModelsBySVtype.TrainRandomForestModelsBySVTYPE as TrainRandomForestModelsBySVTYPE{
                input:
                    training_data = training_data,
                    inheritance = inheritance,
                    genomic_context = genomic_context,
                    svtype = svtype,
                    size_range_list = size_range_list,
                    prefix = prefix,
                    sv_benchmark_docker = sv_benchmark_docker,
                    runtime_attr_override_train_RF_model = runtime_attr_override_train_RF_model
            }
        }


    output{
        Array[Array[File]] pbsv_model = TrainRandomForestModelsBySVTYPE.pbsv_models
        Array[Array[File]] vapor_model = TrainRandomForestModelsBySVTYPE.vapor_models
        
    }
}

