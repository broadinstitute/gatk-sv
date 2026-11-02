# GATK-SV

A structural variation discovery pipeline for Illumina short-read whole-genome sequencing (WGS) data.

For technical documentation on GATK-SV, including how to run the pipeline, please refer to our [website](https://broadinstitute.github.io/gatk-sv/).

## Repository structure
* `/.github`: Continuous integration (CI) and continuous delivery (CD) workflows
* `/dockerfiles`: Resources for building pipeline Docker images
* `/inputs`: Files for generating workflow inputs
  * `/templates`: Input JSON file templates
  * `/values`: Input values used to populate templates
* `/scripts`: Scripts for running tests, building Docker, and analyzing Cromwell metadata files
* `/src`: Main pipeline scripts
  * `/RdTest`: Scripts for depth testing
  * `/sv-pipeline`: Various scripts and packages used throughout the pipeline
  * `/svqc`: Python module for checking that pipeline metrics fall within acceptable limits
  * `/svtest`: Python module for generating various summary metrics from module outputs
  * `/svtk`: Python module of tools for SV-related datafile parsing and analysis
  * `/WGD`: Whole-genome dosage score scripts
* `/wdl`: WDLs running the pipeline. There is a master WDL for running each module, e.g., `ClusterBatch.wdl`.
* `/website`: Website code
