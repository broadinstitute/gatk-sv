---
title: Docker images
description: Docker images
sidebar_position: 4
---

GATK-SV utilizes a set of [Docker](https://www.docker.com/) images for execution within containerized environments.

### Publishing and availability

Dockers are automatically built and pushed to the `us.gcr.io/broad-dsde-methods/gatk-sv` repository under two different conditions:
1. **Release**: upon releasing a new version of GATK-SV. These Dockers are made permanently available.
2. **Commit**: upon merging a new commit to the development branch. These Dockers are ephemeral and may be periodically 
deleted. Any users needing to preserve access to these Docker images should copy them to their own repository. Also
note that these images are built on an "as-needed" basis, meaning an image is only updated if it has been changed
in any way by the commit.

The full set of current images are automatically published and pushed to a public
[Google Artifact Registry](https://cloud.google.com/artifact-registry/docs) and listed in the 
[dockers.json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json)
file.

:::info
Microsoft Azure mirrors of all images can be found in
[dockers_azure.json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers_azure.json),
but these are not currently available for public use.
:::

### Regions (IMPORTANT)

Users using Google Compute VMs outside of the `us-central1` region, i.e. in other NA regions or continents, must copy all
docker images to a repository hosted in their region. If using a Terra workspace, the region is listed under `Cloud 
Information`:`Location` in the workspace dashboard. Please see 
[this article](https://support.terra.bio/hc/en-us/articles/4408985788187-How-to-configure-Google-Artifact-Registry-to-prevent-data-transfer-egress-charges)
for more information.

:::warning
Failure to localize Docker images to your region will incur significant egress costs. 
:::

### Versioning

All Docker images are tagged with a date and version number that must be run with the corresponding version of the 
WDLs. The Docker images built with a particular version can be determined from the `dockers.json` file by checking out
the commit or release of interest and examining `dockers.json`, e.g.
[v0.29-beta](https://github.com/broadinstitute/gatk-sv/blob/v0.29-beta/inputs/values/dockers.json).

Note that a given commit may contain a mixture of Docker image version tags if only a subset of images has actually 
been updated, reflecting the "as-needed" rebuild policy described above.

All Terra workspace releases have the WDLs and Docker images synchronized.

:::info
Cloning a Terra workspace copies a snapshot of the current workspace. Any future updates to the original workspace 
will not be propagated to the clone.
:::

:::warning
We strongly recommend using a single version of GATK-SV for all projects. Switching versions may lead to batch effects 
and/or compromise call set fidelity.
:::


### Workflow inputs

All GATK-SV WDLs utilize a subset of the Docker images, which are specified as inputs. For example, 
the `sv_pipeline` docker is supplied to the input `sv_pipeline_docker` in many WDLs.


### Builds

An in-depth guide to the GATK-SV Docker system and building Docker images is available in the 
[Advanced section](/docs/category/docker-builds).
