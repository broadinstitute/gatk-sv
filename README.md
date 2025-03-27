# GATK-SV

A structural variation discovery pipeline for Illumina short-read whole-genome sequencing (WGS) data.

For technical documentation on GATK-SV, including how to run the pipeline, please refer to our [website](https://broadinstitute.github.io/gatk-sv/).

## Repository structure
* `/dockerfiles`: Resources for building pipeline docker images
* `/inputs`: files for generating workflow inputs
  * `/templates`: Input json file templates
  * `/values`: Input values used to populate templates
* `/scripts`: scripts for running tests, building dockers, and analyzing cromwell metadata files
* `/src`: main pipeline scripts
  * `/RdTest`: scripts for depth testing
  * `/sv-pipeline`: various scripts and packages used throughout the pipeline
  * `/svqc`: Python module for checking that pipeline metrics fall within acceptable limits
  * `/svtest`: Python module for generating various summary metrics from module outputs
  * `/svtk`: Python module of tools for SV-related datafile parsing and analysis
  * `/WGD`: whole-genome dosage score scripts
* `/wdl`: WDLs running the pipeline. There is a master WDL for running each module, e.g. `ClusterBatch.wdl`.
* `/website`: website code

## CI/CD
This repository is maintained following the norms of 
continuous integration (CI) and continuous delivery (CD). 
GATK-SV CI/CD is developed as a set of Github Actions
workflows that are available under the `.github/workflows`
directory. Please refer to the [workflow's README](.github/workflows/README.md) 
for their current coverage and setup. 