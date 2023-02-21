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

# Cromwell

[Cromwell](https://github.com/broadinstitute/cromwell) is a workflow management system
that takes a workflow (e.g., a workflow written in [Workflow Description Language (WDL)](https://openwdl.org)), 
its dependencies and input data, and runs it on a given platform 
(e.g., [GCP](https://cromwell.readthedocs.io/en/stable/backends/Google/)). 
In order to run a workflow on Cromwell, you need a running instance of 
Cromwell that is available in two forms: [Server and stand-alone mode](https://cromwell.readthedocs.io/en/stable/Modes/).

In general, you may use a managed Cromwell server maintained by your 
institute or host a self-managed server, or run Cromwell as a standalone Java application.
The former is ideal for large scale execution in a managed environment, 
and the latter is useful for small scale and isolated WDL development.

:::info
Due to its dependency on cloud-hosted resources and large-scale execution needs,
we currently do not support running the entire GATK-SV pipeline using 
Cromwell as a [stand-alone Java application](https://cromwell.readthedocs.io/en/stable/Modes/#run) 
:::

# Cromwell Server

There are two option to communicate with a running Cromwell server: 
[REST API](https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/), and
[Cromshell](https://github.com/broadinstitute/cromshell) which is a command line tool
to interface with a Cromwell server. We recommend using Cromshell due to its simplicity 
of use. This documentation is explained using Cromshell, but the same steps can be 
taken using the REST API.

- **Setup Cromwell**: You may follow [this](https://cromwell.readthedocs.io/en/stable/Modes/) documentation
  on setting up a Cromwell server.

- **Setup Cromshell**: You may follow [this](https://github.com/broadinstitute/cromshell) documentation
  on installing and configuring Cromshell to communicate with the Cromwell server. 
