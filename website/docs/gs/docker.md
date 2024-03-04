---
title: Docker Images
description: GATK-SV Docker Images
sidebar_position: 4
slug: ./dockers
---


To make the analysis process scalable, reproducible, and cost-efficient,
GATK-SV is designed as a cloud-native pipeline, 
meaning it runs on virtual machines (VMs) in the cloud, 
which are pre-configured with all the necessary tools, scripts, 
and settings for reliable analysis. To easily replicate and share 
the analysis, GATK-SV uses Docker technology. Docker packages the tools, 
scripts, and their requirements into self-contained units called containers. 
These containers can be deployed on different VMs in the cloud, 
ensuring consistent and reproducible analysis for various experiments 
and collaborations.

The latest Docker image builds can be found in the following files.



- [`dockers.json`](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json).
  The list of images hosted on Google Container Registry (GCR). 
  You may use the Docker images listed in this file if you are running 
  the pipeline on Google Cloud Platform (GCP).

- [`dockers_azure.json`](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers_azure.json).
  The list of images hosted on Azure Container Registry (ACR).
  You may use the Docker images listed in this file if you are
  running the pipeline on Azure.


:::tip For developers and power users

You may refer to [this section](/docs/advanced/docker/) for a detailed 
description of the Docker images, including their design principles, 
as well as guides on build and deploy them.
:::
  