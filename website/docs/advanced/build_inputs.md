---
title: Building inputs
description: Building work input json files
sidebar_position: 1
slug: build_inputs
---

Each workflow of the GATK-SV pipeline takes a unique set of arguments as inputs. 
You have different options for configuring them depending on the platform you're 
using to run the pipeline. 
For instance, you may use Terra workspaces if you're running on Terra (user-friendly), 
or JSON files if you're running on Cromwell (for development and advanced use-cases). 
For each workflow, we provide example configurations that help both in setting up 
your own Terra workspace or for testing purposes with sample data. 
You may run the following commands to get these example inputs.


1. Clone GATK-SV (you may skip this step if you have already done so).

    ```shell
   git clone https://github.com/broadinstitute/gatk-sv && cd gatk-sv
    ```

2. Create a JSON file containing the Terra billing project (for use on Terra) 
   or the Google project ID (for use on Cromwell) that you will use to run 
   the workflows with the test input. You may create this file by running
   the following command and replacing `"my-google-project-id"` and 
   `"my-terra-billing-project"` with your project and billing IDs. 

    ```shell
    echo '{ "google_project_id": "my-google-project-id", "terra_billing_project_id": "my-terra-billing-project" }' > inputs/values/google_cloud.my_project.json
    ```

3. Create test inputs.

    ```shell
    bash scripts/inputs/build_default_inputs.sh -d . -c google_cloud.my_project 
    ```
   
    Running this command generates test inputs in `gatk-sv/inputs/build` with the following structure. 

   ```shell
   inputs/build
   ├── NA12878
   │   ├── terra
   │   └── test
   ├── NA19240
   │   └── test
   ├── hgdp
   │   └── test
   └── ref_panel_1kg
       ├── terra
       └── test
   ```

## Building inputs for specific use-cases (Advanced)

### Build for batched workflows

```shell
python scripts/inputs/build_inputs.py \
  inputs/values \
  inputs/templates/test/GATKSVPipelineSingleSample \
  inputs/build/NA19240/test \
  -a '{ "test_batch" : "ref_panel_1kg", "cloud_env": "google_cloud.my_project" }'
```


### Generating a reference panel

This section only applies to the single-sample mode. 
New reference panels can be generated from a single run of the 
`GATKSVPipelineBatch` workflow. 
If using a Cromwell server, we recommend copying the outputs to a 
permanent location by adding the following option to the 
[workflow configuration](https://cromwell.readthedocs.io/en/latest/wf_options/Overview/) 
file:

```json
"final_workflow_outputs_dir" : "gs://my-outputs-bucket",
"use_relative_output_paths": false,
```

Here is an example of how to generate workflow input jsons from `GATKSVPipelineBatch` workflow metadata:

1. Get metadata from Cromwshell.

   ```shell
   cromshell -t60 metadata 38c65ca4-2a07-4805-86b6-214696075fef > metadata.json
   ```

2. Run the script. 

   ```shell
   python scripts/inputs/create_test_batch.py \
      --execution-bucket gs://my-exec-bucket \
      --final-workflow-outputs-dir gs://my-outputs-bucket \
      metadata.json \
      > inputs/values/my_ref_panel.json
   ```

3. Define your google project id (for Cromwell inputs) and Terra billing project (for workspace inputs).

   ```shell
   echo '{ "google_project_id": "my-google-project-id", "terra_billing_project_id": "my-terra-billing-project" }' > inputs/values/google_cloud.my_project.json
   ```
   
4. Build test files for batched workflows (google cloud project id required).

   ```shell
   python scripts/inputs/build_inputs.py \
      inputs/values \
      inputs/templates/test \
      inputs/build/my_ref_panel/test \
      -a '{ "test_batch" : "ref_panel_1kg", "cloud_env": "google_cloud.my_project" }'
   ```

5. Build test files for the single-sample workflow

   ```shell
   python scripts/inputs/build_inputs.py \
       inputs/values \
       inputs/templates/test/GATKSVPipelineSingleSample \
       inputs/build/NA19240/test_my_ref_panel \
       -a '{ "single_sample" : "test_single_sample_NA19240", "ref_panel" : "my_ref_panel" }'
   ```
   
6. Build files for a Terra workspace.

   ```shell 
   python scripts/inputs/build_inputs.py \
      inputs/values \
      inputs/templates/terra_workspaces/single_sample \
      inputs/build/NA12878/terra_my_ref_panel \
      -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "my_ref_panel" }'
   ```
   
Note that the inputs to `GATKSVPipelineBatch` may be used as resources 
for the reference panel and therefore should also be in a permanent location.
