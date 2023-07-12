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


- **Modular and focused structure:** 
Each image includes task-specific tools, simplifying the use and maintenance of 
GATK-SV Docker images for users and developers, respectively.


- **Reduced Docker image size:**
Using task-specific Docker images reduces sizes, requiring less storage space 
in container registries such as Azure Container Registry (ACR) or 
Google Cloud Container Registry (GCR). It also enables faster image transfer 
when creating virtual machines for task execution.


- **Enhanced maintenance and extensibility:**
Maintainers can easily modify specific tools or configurations within 
a single image without affecting others, improving maintainability and 
facilitating seamless expansion by adding or replacing tools as required.


- **Consistency and efficiency:**
Building images on top of existing setups and tools promotes code 
reuse and reduces duplication, ensuring consistent configurations 
across pipeline stages. It simplifies dependency management by 
allowing changes or updates at the appropriate level, cascading 
down to dependent images.


In summary, splitting tools into smaller, task-specific 
Docker images optimizes storage, execution, maintenance, and extensibility. 
It enhances consistency, code reuse, and dependency management, 
ensuring efficient and scalable pipeline execution.
