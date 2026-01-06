---
title: Building inputs
description: Building work input json files
sidebar_position: 3
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

2. Create test inputs.

    ```shell
    bash scripts/inputs/build_default_inputs.sh -d .
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

### Build for batched workflows

```shell
python scripts/inputs/build_inputs.py \
  inputs/values \
  inputs/templates/test/GATKSVPipelineSingleSample \
  inputs/build/NA19240/test \
  -a '{ "test_batch" : "ref_panel_1kg" }'
```

