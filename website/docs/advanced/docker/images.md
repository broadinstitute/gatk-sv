---
title: Docker Images Hierarchy
description: Docker Image Dependencies
sidebar_position: 1
---

import useBaseUrl from '@docusaurus/useBaseUrl';
import ThemedImage from '@theme/ThemedImage';

:::info
This page provides a detailed explanation of Docker 
images and their hierarchy. For information on the process 
of building these images, please refer to [this section](/docs/advanced/docker/deploy).
:::


The tools, scripts, dependencies, and configurations utilized by the 
GATK-SV pipeline, written in WDL, are organized into separate Docker 
containers. This modular approach ensures that each container 
contains only the necessary tools for its specific task, 
resulting in smaller image sizes. This design choice simplifies 
the definition of Dockerfiles and facilitates easier maintenance. 
Moreover, the smaller image sizes contribute to reduced disk 
usage and lower workflow execution costs.


The figure below illustrates the relationships between the GATK-SV Docker images.


<ThemedImage
  alt="Interconnection Between GATK-SV Docker Images"
  sources={{
    light: useBaseUrl('/img/docker_hierarchy.png'),
    dark: useBaseUrl('/img/docker_hierarchy.png'),
  }}
/>

The image depicts the hierarchical relationship among GATK-SV 
Docker images. Arrows indicate the flow from a base image 
to a derived image. The base image, located at the arrow's 
starting point, shares its content which is then expanded 
upon and modified in the derived image. In simple terms, 
the derived image inherits the same tools and configuration 
as the base image, while incorporating additional settings and tools.


The list of the Docker images and their latest builds 
are available in [`dockers.json`](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json)
and [`dockers_azure.json`](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers_azure.json)
for images hosted on Google Container Registry (GCR) and Azure Container Registry (ACR), respectively.


## Advantages of Dividing Images by Functionality

The GATK-SV pipeline utilizes Docker images to encapsulate the necessary tools, 
dependencies, and configurations. Instead of having a single monolithic image, 
the pipeline is organized into multiple smaller images, each focusing on a specific task. 
This approach offers several benefits.


By splitting the tools into separate Docker images, we achieve a modular 
and focused structure. Each image contains the tools required for a specific 
task within the GATK-SV pipeline. This enables users and developers to easily
work with individual images, as they can identify the specific tools needed 
for their particular analysis.


Moreover, using smaller, task-specific Docker images offers the advantage 
of reduced sizes, which is particularly beneficial in cloud environments. 
These smaller images require less storage space when stored in container 
registries like Google Cloud Container Registry (GCR) or Azure Container Registry (ACR). 
Additionally, when creating virtual machines for workflow task execution, 
the transfer of these smaller images is more efficient.


Separate Docker images enhance maintenance and extensibility 
in the GATK-SV pipeline. Maintainers can easily modify or update 
specific tools or configurations within a single image without 
impacting others. This granularity improves maintainability 
and enables seamless expansion of the pipeline by adding or 
replacing tools as required.


Additionally, the Docker image hierarchy offers advantages in terms of 
consistency and efficiency. One image can be built upon another, 
leveraging existing setups and tools. This promotes code reuse and 
reduces duplication, resulting in consistent configurations across 
different stages of the pipeline. It also simplifies the management 
of common dependencies, as changes or updates can be applied at the 
appropriate level, cascading down to the dependent images.


In summary, by splitting the tools into smaller, task-specific images, 
the pipeline becomes more modular and manageable. 
This approach optimizes storage, execution, maintenance, 
and extensibility in cloud environments. 
Leveraging Docker's image hierarchy further enhances consistency, 
code reuse, and dependency management, ensuring efficient and 
scalable execution of the pipeline.
