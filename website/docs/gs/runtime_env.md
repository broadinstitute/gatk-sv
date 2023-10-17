---
title: Runtime Environments
description: Describes the supported runtime environments.
sidebar_position: 7
slug: ./runtime-env
---

The GATK-SV pipeline consists of _workflows_ and _reference data_ that
orchestrates the analysis flow of input data. Hence, a successful
execution requires running the _workflows_ on _reference_ and input data.

:::info Currently supported backends: GCP
GATK-SV has been tested only on the Google Cloud Platform (GCP); 
therefore, we are unable to provide specific guidance or support 
for other execution platforms including HPC clusters and AWS.
:::

## Alternative backends

Contributions from the community to improve portability between backends 
will be considered on a case-by-case-basis. We ask contributors to 
please adhere to the following guidelines when submitting issues and pull requests:

1. Code changes must be functionally equivalent on GCP backends, i.e. not result in changed output
2. Increases to cost and runtime on GCP backends should be minimal
3. Avoid adding new inputs and tasks to workflows. Simpler changes 
   are more likely to be approved, e.g. small in-line changes to scripts or WDL task command sections
4. Avoid introducing new code paths, e.g. conditional statements
5. Additional backend-specific scripts, workflows, tests, and Dockerfiles will not be approved
6. Changes to Dockerfiles may require extensive testing before approval

We still encourage members of the community to adapt GATK-SV for non-GCP backends 
and share code on forked repositories. Here are a some considerations:

- Refer to Cromwell's [documentation](https://cromwell.readthedocs.io/en/stable/backends/Backends/) 
  for configuration instructions.

- The handling and ordering of `glob` commands may differ between platforms.

- Shell commands that are potentially destructive to input files 
  (e.g. `rm`, `mv`, `tabix`) can cause unexpected behavior on shared filesystems. 
  Enabling [copy localization](https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options) 
  may help to more closely replicate the behavior on GCP.

- For clusters that do not support Docker, Singularity is an alternative. 
  See [Cromwell documentation on Singularity](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#singularity).

- The GATK-SV pipeline takes advantage of the massive parallelization possible in the cloud. 
  Local backends may not have the resources to execute all of the workflows. 
  Workflows that use fewer resources or that are less parallelized may be more successful. 
  For instance, some users have been able to run [GatherSampleEvidence](#gather-sample-evidence) on a SLURM cluster.
