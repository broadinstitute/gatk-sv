---
title: Runtime environments
description: Describes the supported runtime environments.
sidebar_position: 5
slug: ./runtime-env
---

The GATK-SV pipeline consists of workflows implemented in the Workflow Description Language 
([WDL](https://openwdl.org/)) and is built for use on the Google Cloud Platform ([GCP](https://cloud.google.com/)). 

### Terra (recommended)
To facilitate easy-of-use, security, and collaboration, GATK-SV is available on the [Terra](https://app.terra.bio/) 
platform. Users should clone pre-configured Terra workspaces as a starting point for running GATK-SV:

- [Single-sample workspace](https://app.terra.bio/#workspaces/help-gatk/GATK-Structural-Variants-Single-Sample)
- [Joint calling workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling)

These workspaces are actively maintained by the development team and will be updated with critical fixes and major releases.

### Cromwell (advanced)
Advanced users and developers who wish to run GATK-SV on a dedicated Cromwell server using GCP as a backend should refer 
to the [Advanced Guides](/docs/category/advanced-guides) section.

### Alternative backends (not supported)

Use of other backends (e.g. AWS or on-prem HPC clusters) is not currently supported. However, contributions from the 
community to improve portability between backends will be considered on a case-by-case-basis. We ask contributors to 
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
  For instance, some users have been able to run [GatherSampleEvidence](../modules/gse) on a SLURM cluster.
