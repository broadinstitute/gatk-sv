---
title: Overview
description: Introduction to Cromwell
sidebar_position: 0
---

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
Cromwell as a [stand-alone Java application](https://cromwell.readthedocs.io/en/stable/Modes/#run).

Additionally, we currently only support running the pipeline on 
Google Cloud Platform (GCP).
:::

# Cromwell Server

There are two options to communicate with a running Cromwell server: 
[REST API](https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/), and
[Cromshell](https://github.com/broadinstitute/cromshell) which is a command line tool
to interface with a Cromwell server. We recommend using Cromshell due to its simplicity 
of use. This documentation is explained using Cromshell, but the same steps can be 
taken using the REST API.

- **Setup Cromwell**: You may follow [this](https://cromwell.readthedocs.io/en/stable/Modes/) documentation
  on setting up a Cromwell server.

- **Setup Cromshell**: You may follow [this](https://github.com/broadinstitute/cromshell) documentation
  on installing and configuring Cromshell to communicate with the Cromwell server. 