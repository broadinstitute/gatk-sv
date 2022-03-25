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

The output of this step is a single VCF with a `BS` attribute added to the sample-level genotype fields for all variants.    

## Step 3: Filter VCF based on lgBoost scores  

This step filters low-quality genotypes based on their Boost scores as annotated in `Step 2`.  

Use `wdl/Module07FilterGTsPart3ApplyModel.wdl` to run this step.  

This step requires a minimum BS lookup table for filtering. This table can be generated in several ways, such as:  

1. Manually generating a table of hard thresholds; see `example_inputs/flat_BS_filter_lookup_table.tsv` for an example of this format, or
2. Dynamically optimizing thresholds for different types of SVs using the minGQ workflows. Note that this option is still in development.  

The output of this step is a filtered VCF ready for additional downstream filtering & annotation.  

---  

## OLD (pre-RLC):  

## Layer 2: Integrate boost score across samples and calculate the overall boost failure rate for each SV locus.
This step integrates the boost score assigned to each sample in Layer 1. A boost score cutoff should be decided / arbitrarily assigned for each variant type. For a specific SV, samples that carry this SV but have a lower boost score than the cutoff would be considered as failure, and vise versa.

Use `IntegrateBoostScores.wdl` to run this step. 

Example of the json: *./example_jsons/IntegrateBoostResultsAcrossBatches.bscff_minus_1.json*

## others
Example of the boost score cutoff table can be found at: *./example_boost_cff_table/boost_cff_table_minus1.tsv*

The docker files for `sv_benchmark_docker` can be found at: gatk-sv/dockerfiles/rdpesr_docker

