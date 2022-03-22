## _Post hoc_ genotype filtering with lgBoost

This directory includes scripts to train and apply a light gradient-boosted (lgBoost) model for genotype filtering.  

A trained lgBoost model can be applied with the following steps:  

### Step 1: Apply lgBoost to each sample  

This step applies the trained lgBoost model on each individual genome and assigns a boost score for each SV per sample.  

Higher boost score indicates higher genotype quality, and vice versa.

Use `TrainAndApplyBoostFilter.wdl` to run this step. 

Example json: `./example_jsons/TrainAndApplyBoostModel.group_all.json`  

This step splits all samples into batches, applies the boost model to all variants per sample, and returns one tarball per batch.  

Each tarball contains a two-column .tsv per sample with the variant ID and boost score for every non-reference variant per sample.

## Step 2: Annotate VCF with lgBoost results  

This step aggregates the results from `Step 1` and writes them into the input VCF under a custom genotype field, `BS`.  

Use `AnnotateVcfWithBoostScores.wdl` to run this step.  

This step is currently a work-in-progress. Documentation will be updated soon.  

---  

## OLD (pre-RLC):  

## Layer 2: Integrate boost score across samples and calculate the overall boost failure rate for each SV locus.
This step integrates the boost score assigned to each sample in Layer 1. A boost score cutoff should be decided / arbitrarily assigned for each variant type. For a specific SV, samples that carry this SV but have a lower boost score than the cutoff would be considered as failure, and vise versa.

Use `IntegrateBoostScores.wdl` to run this step. 

Example of the json: *./example_jsons/IntegrateBoostResultsAcrossBatches.bscff_minus_1.json*

## others
Example of the boost score cutoff table can be found at: *./example_boost_cff_table/boost_cff_table_minus1.tsv*

The docker files for `sv_benchmark_docker` can be found at: gatk-sv/dockerfiles/rdpesr_docker

