---
title: Quick Start
description: Run the pipeline on demo data.
sidebar_position: 1
slug: ./qs
---

This page provides steps for running the pipeline using demo data. 

# Quick Start on Terra

https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/gatk-sv-pipeline-1kgp

# Quick Start on Cromwell

This section walks you through the steps of running pipeline using 
demo data on a managed Cromwell server.

### 1/4: Setup Environment

- A running instance of a Cromwell server. While you can 
  [host your own Cromwell server](https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/),
  (or [these steps](cromwell#setup-cromwell)) due to the required domain expertise, we recommend using a managed 
  Cromwell instance for the purpose of _Quick Start_.

- Install Cromshell and configure it to connect with the Cromwell server you are using. 
  You may follow [these steps](cromwell#setup-cromshell).

### 2/4: Build Inputs

- Example workflow inputs can be found in `/inputs`. 
  Build using `scripts/inputs/build_default_inputs.sh`, 
  which generates input jsons in `/inputs/build`.

- Some workflows require a Google Cloud Project ID to be defined in 
  a cloud environment parameter group. Workspace builds require a 
  Terra billing project ID as well. An example is provided at 
  /inputs/values/google_cloud.json but should not be used, 
  as modifying this file will cause tracked changes in the repository. 
  Instead, create a copy in the same directory with the format 
  google_cloud.my_project.json and modify as necessary.

  Note that these inputs are required only when certain data are 
  located in requester pays buckets. If this does not apply, 
  users may use placeholder values for the cloud configuration 
  and simply delete the inputs manually.

### 3/4: Requester Pays Buckets

The following parameters must be set when certain input data is in requester pays (RP) buckets:

`GATKSVPipelineSingleSample.requester_pays_cram` and 
`GATKSVPipelineBatch.GatherSampleEvidenceBatch.requester_pays_crams` - 
set to `True` if inputs are CRAM format and in an RP bucket, otherwise `False`.


### 4/4: Execution

```shell
> mkdir gatksv_run && cd gatksv_run
> mkdir wdl && cd wdl
> cp $GATK_SV_ROOT/wdl/*.wdl .
> zip dep.zip *.wdl
> cd ..
> echo '{ "google_project_id": "my-google-project-id", "terra_billing_project_id": "my-terra-billing-project" }' > inputs/values/google_cloud.my_project.json
> bash scripts/inputs/build_default_inputs.sh -d $GATK_SV_ROOT -c google_cloud.my_project
> cp $GATK_SV_ROOT/inputs/build/ref_panel_1kg/test/GATKSVPipelineBatch/GATKSVPipelineBatch.json GATKSVPipelineBatch.my_run.json
> cromshell submit wdl/GATKSVPipelineBatch.wdl GATKSVPipelineBatch.my_run.json cromwell_config.json wdl/dep.zip
```

where `cromwell_config.json` is a Cromwell 
[workflow options file](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/).
Note users will need to re-populate batch/sample-specific parameters (e.g. BAMs and sample IDs).
