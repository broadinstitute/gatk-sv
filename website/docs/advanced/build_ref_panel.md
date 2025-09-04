---
title: Building reference panels
description: Building reference panels for the single-sample pipeline
sidebar_position: 4
slug: build_ref_panel
---

## Overview

The [single-sample mode](/docs/gs/calling_modes#single-sample-mode) requires a reference panel of samples. GATK-SV provides a default reference panel described [here](/docs/execution/single#reference-panel). However, users may wish to generate new reference panels. The panel samples should be chosen using the same criteria as [joint calling batching](/docs/modules/eqc#batching), i.e. even sex balance and similar depth, WGD score, insert size, aligner, etc. with respect to each other and with respect to all samples that will be called in single-sample mode. The reference panel may only consist of a single batch - multi-batch reference panels are not supported. 

In addition, we recommend the following guidelines:

- Maximize genetic diversity
- No related samples
- No aneuploidies, germline or mosaic, on allosomes or autosomes
- No samples with low sequencing coverage
- No samples in tails of WGD score distribution

### Run joint calling

First run the new reference panel samples through the [joint calling Terra workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling). All numbered workflows must be run. Refer to the [joint calling documentation](/docs/execution/joint) for further instructions. Optional workflows without numbers do not need to be run.

:::important
The default configurations should be used for all workflows. The following steps assume that specific output fields are populated in the workspace data table and that workspace attributes are set according to the default configuration. However, users may adjust numeric parameters such as genotype filtering cutoffs.
:::

### Run notebook

The joint calling workspace contains a Jupyter notebook, `CreateReferencePanel.ipynb`, that should next be run to generate json-encoded resources that are required for building inputs for single-sample calling. Navigate to this notebook in the `Analyses` section of the joint calling workspace and follow the instructions.

### Clone single-sample Terra workspace and check version

Before continuing, you must determine which version of the single-sample pipeline you will be running. To do that, create a clone of the [single-sample Terra workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Single-Sample). 

Next, inspect the box containing `gatk-sv-single-sample` to find the current version after the `V`. For example:

```
gatk-sv-single-sample

V v0.26.9-beta
Source: Dockstore
```
indicates that the current version is `v0.26.9-beta`. Take note of the version, as you will need it in the next step.

:::warning
The public Terra workspaces are kept up to date with the latest versions of joint calling and single-sample modes that have undergone testing. The versions may be out of sync, however, but should not generally be mixed, and the above step ensures that the latest safe version is being used. Users may also elect newer versions through the workflow configuration, but be aware that it may not be fully tested.
:::

### Clone and checkout Git repository

The notebook will have generated a json file that can be consumed by the GATK-SV [inputs generation framework](/docs/advanced/build_inputs). Create a clone of the Git repository and checkout the current version:

```shell
github clone https://github.com/broadinstitute/gatk-sv.git
cd gatk-sv
git checkout RELEASE_VERSION
```

where `RELEASE_VERSION` is the workflow version from the cloned workspace in the previous step.

### Download resources json

The path to the reference panel resources json is printed out after the last cell, for example:

```
File test_panel.json uploaded to gs://fc-7dd8986b-d916-46b0-ba1a-8b09f80f7b83/json/test_panel.json
```

Download the json to the `inputs/values/` subdirectory in your gatk-sv git clone:

```shell
gsutil cp gs://fc-7dd8986b-d916-46b0-ba1a-8b09f80f7b83/json/test_panel.json ./inputs/values/
```

### Build single-sample Terra workspace configuration

Next run the following command from the root of your local gatk-sv Git clone:

```shell 
python scripts/inputs/build_inputs.py \
  inputs/values \
  inputs/templates/terra_workspaces/single_sample \
  inputs/build/NA12878/MY_TERRA_CONFIG \
  -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "REF_PANEL_NAME" }'
```

where `REF_PANEL_NAME` again exactly matches the corresponding variable from the notebook (the json is named `REF_PANEL_NAME.json`). In addition, `MY_TERRA_CONFIG` can be renamed if desired.

Confirm that the build was successful by running:

```shell
ls inputs/build/NA12878/MY_TERRA_CONFIG
```

and seeing the output:

```
GATKSVPipelineSingleSample.json
participant.tsv
sample.tsv
single_sample_workspace_dashboard.md
workspace.tsv
```

If the directory does not exist then the build was not successful. If this occurs, run the `build_inputs.py` script with the `--log-info` flag to print troubleshooting logs.

### Create and configure Terra workspace {#configure-terra}

Now we will configure a Terra workspace with the reference panel in the following steps:

1. Return to your clone of the single-sample Terra workspace. 
2. Navigate to the `Data` tab and click `Workspace Data` on the left navigation bar.
3. It is good practice to clear all existing entries. To do this, click on the top-left checkbox to select all rows, then click `Edit` and `Delete selected variables`. 
4. Populate the workspace attributes with the `workspace.tsv` file built in the previous section, either through the `Import Data` wizard or by dragging and dropping the file into your browser.
5. Navigate to the `gatk-sv-single-sample` workflow configuration in the `Workflows` tab. 
6. Reset the current configuration by clicking `Clear inputs`.
7. Update the inputs with the `GATKSVPipelineSingleSample.json` file built in the last step, either by clicking `upload json` or dragging and dropping it in to your browser. 
8. Click `Save` to commit the update.

