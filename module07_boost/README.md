This repository includes scripts to train and apply boost model for variant filtering.

The lgBoost model should be applied in layers, below are the instructions.

## Layer 1: Filter SVs in each individual genome.
This step applies the trained lgBoost model on each individual genome and assigns a boost score for each SV in each genome. Higher boost score indicates better quality of the variants, and vise versa.

Use `TrainAndApplyBoostModel.wdl` to run this step. Example of the json: *./example_jsons/TrainAndApplyBoostModel.group_all.json*

## Layer 2: Integrate boost score across samples and calculate the overall boost failure rate for each SV locus.
This step integrates the boost score assigned to each sample in Layer 1. A boost score cutoff should be decided / arbitrarily assigned for each variant type. For a specific SV, samples that carry this SV but have a lower boost score than the cutoff would be considered as failure, and vise versa.

Use `IntegrateBoostScores.wdl` to run this step. 

Example of the json: *./example_jsons/IntegrateBoostResultsAcrossBatches.bscff_minus_1.json*

Example of the boost score cutoff table can be found at: *./example_boost_cff_table/boost_cff_table_minus1.tsv*

The docker files for `sv_benchmark_docker` can be found at: gatk-sv/dockerfiles/rdpesr_docker

