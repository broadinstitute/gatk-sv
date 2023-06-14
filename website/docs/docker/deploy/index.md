---
title: Deploying Docker Images
description: Docker Concepts and Execution Overview
sidebar_position: 2
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

:::info
This section offers a comprehensive explanation of the process of 
building, testing, and publishing Docker images. For details 
regarding the images and their hierarchy, please refer to
[this page](/docs/docker/images).
:::


GATK-SV Docker image _deployment_ involves the essential steps of 
_building_, _testing_, and _publishing_ to Docker container registries. 
There are two deployment options available: fully automated and manual. 
With the fully automated approach, GATK-SV Docker images are built 
and published to Google Container Registry (GCR) and 
Azure Container Registry (ACR) through continuous integration and 
continuous delivery (CI/CD) after merging a pull request. 
However, if you are working on extending or improving the 
GATK-SV Docker images, you may need to build the images locally 
for testing or store them on an alternative container registry. 
This section provides comprehensive insights into the automatic 
build process and a detailed guide on locally building the images 
for development purposes.
