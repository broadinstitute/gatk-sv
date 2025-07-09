---
title: Building reference panels
description: Building reference panels for the single-sample pipeline
sidebar_position: 4
slug: build_ref_panel
---

## Overview

The [single-sample mode](/docs/gs/calling_modes#single-sample-mode) requires a reference panel of samples. GATK-SV provides a default reference panel described [here](/docs/execution/single#reference-panel). However, users may wish to generate new reference panels. The panel samples should be chosen using the same criteria as [joint calling batching](/docs/modules/eqc#batching), i.e. even sex balance and similar depth, WGD score, insert size, aligner, etc. with respect to each other and with respect to all samples that will be called in single-sample mode.

### Run joint calling

First run the new reference panel samples through the [joint calling Terra workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling). All numbered workflows must be run. Refer to the [joint calling documentation](/docs/execution/joint) for further instructions. Optional workflows without numbers do not need to be run.

:::important
The default configurations should be used for all workflows. The following steps assume that specific output fields are populated in the workspace data table and that workspace attributes are set according to the default configuration.
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
The public Terra workspace is kept up to date with the latest version that has undergone testing. Users may elect a newer version through the workflow configuration, but be aware that it may not be fully tested.
:::

### Clone and checkout Git repository

The notebook will have generated a json file that can be consumed by the GATK-SV [inputs generation framework](/docs/advanced/build_inputs). Create a clone of the Git repository and checkout the current version:

```shell
> github clone https://github.com/broadinstitute/gatk-sv.git
> cd gatk-sv
> git checkout RELEASE_VERSION
```

where `RELEASE_VERSION` is the workflow version from the cloned workspace in the previous step.

### Create resources json

Next create an empty file to populate with the reference panel resources:

```shell
> touch ./inputs/values/REF_PANEL_NAME.json
```

The file name should exactly match the `REF_PANEL_NAME` variable from the notebook. Open the new json in a text editor and copy and paste the notebook output into the file. Save the file and exit your editor.

### Build single-sample Terra workspace configuration

Next run the following command from the root of your local gatk-sv Git clone:

   ```shell 
   > python scripts/inputs/build_inputs.py \
      inputs/values \
      inputs/templates/terra_workspaces/single_sample \
      inputs/build/NA12878/MY_TERRA_CONFIG \
      -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "REF_PANEL_NAME" }'
   ```

where `REF_PANEL_NAME` again exactly matches the corresponding variable from the notebook. In addition, `MY_TERRA_CONFIG` can be renamed if desired.

Confirm that the build was successful:

   ```shell
   > ls inputs/build/NA12878/MY_TERRA_CONFIG
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

