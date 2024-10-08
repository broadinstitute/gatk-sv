---
title: Run
description: Running GATK-SV on Cromwell
sidebar_position: 1
slug: ./qs
---

# Quick Start on Cromwell

This section walks you through the steps of running pipeline using 
demo data on a managed Cromwell server.

### Environment Setup

- A running instance of a Cromwell server.

- Install Cromshell and configure it to connect with the Cromwell server you are using. 
  You may refer to the documentation on [Cromshell README](https://github.com/broadinstitute/cromshell).

### Build Inputs

We provide options for building example inputs that you may use as a reference 
to configure a Terra workspace or Cromwell submissions (advanced) with your own data. 
Please refer to [this page](/docs/advanced/build_inputs) for instructions on how to build these inputs.

### MELT
Important: The example input files contain MELT inputs that are NOT public 
(see [Requirements](https://github.com/broadinstitute/gatk-sv#requirements)). These include:

- `GATKSVPipelineSingleSample.melt_docker` and `GATKSVPipelineBatch.melt_docker` - MELT docker URI 
(see [Docker readme](https://github.com/talkowski-lab/gatk-sv-v1/blob/master/dockerfiles/README.md))
- `GATKSVPipelineSingleSample.ref_std_melt_vcfs` - Standardized MELT VCFs ([GatherBatchEvidence](/docs/modules/gbe))
The input values are provided only as an example and are not publicly accessible. 
- In order to include MELT, these values must be provided by the user. MELT can be 
  disabled by deleting these inputs and setting `GATKSVPipelineBatch.use_melt` to false.

### Execution

```shell
> mkdir gatksv_run && cd gatksv_run
> mkdir wdl && cd wdl
> cp $GATK_SV_ROOT/wdl/*.wdl .
> zip dep.zip *.wdl
> cd ..
> bash scripts/inputs/build_default_inputs.sh -d $GATK_SV_ROOT
> cp $GATK_SV_ROOT/inputs/build/ref_panel_1kg/test/GATKSVPipelineBatch/GATKSVPipelineBatch.json GATKSVPipelineBatch.my_run.json
> cromshell submit wdl/GATKSVPipelineBatch.wdl GATKSVPipelineBatch.my_run.json cromwell_config.json wdl/dep.zip
```

where `cromwell_config.json` is a Cromwell 
[workflow options file](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/).
Note users will need to re-populate batch/sample-specific parameters (e.g. BAMs and sample IDs).
