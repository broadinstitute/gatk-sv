---
title: Manual Deployment
description: Build and Publish Images
sidebar_position: 3
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';


If you are contributing to the GATK-SV codebase, specifically focusing on 
enhancing tools, configuring dependencies in Dockerfiles, or modifying GATK-SV scripts 
within the Docker images, it is important to build and test the Docker images locally. 
This ensures that the images are successfully built and function as intended.

The process of updating GATK-SV Docker images involves two steps: build and publish.

- **Build**: Create Docker images from Dockerfiles and store them on your computer.

- **Publish**: Upload the built Docker images to container registries 
(e.g., Google Container registry, or Azure container registry) 
to make them available for use in Terra or Cromwell.

You may refer to [this page](/docs/advanced/docker/index.md) for detailed description of the process. 
To streamline the process, we have developed a Python script 
that automates the image building and publishing to your container registry.
This section provides guidelines on building and publishing the images using this script. 


:::warning Linux Machine Required

Only Linux machines (dedicated or virtual) are supported for building GATK-SV Docker images. 
Images created on non-Linux machines may not work with Terra or Cromwell execution environment.
The instructions provided on this page assume you are using a Linux Ubuntu machine.
:::



## Setup

### Runtime environment {#runtime}

Currently, GATK-SV Docker images can only be built on the `linux/amd64` platform, 
which is a machine running Linux OS on x86-64 architecture.
Images build on Apple M1 (`linux/arm64`) are not currently supported.
You can use a local Linux machine or obtain a virtual machine from a cloud platform. 

You may follow the steps in the 
[GCP](https://cloud.google.com/compute/docs/instances/create-start-instance#publicimage)
or [Azure](https://learn.microsoft.com/en-us/azure/virtual-machines/windows/quick-create-portal)
documentation to create a virtual machine (VM) on Google Cloud Platform (GCP) or Microsoft Azure respectively. 
Make sure the VM is built using an Ubuntu image, has at least 8 GB RAM, and some additional 
disk space (e.g., 50 GB should be sufficient).

Building and publishing GATK-SV Docker images is time-consuming and can take around 1 hour.
Therefore, we recommend using a terminal multiplexer 
(e.g., [tmux](https://github.com/tmux/tmux/wiki/Getting-Started); 
[tmux cheat sheet](https://tmuxcheatsheet.com)) 
when running on a VM to ensure the process continues even if you are disconnected from the VM.

### Docker {#docker}

[Install](https://docs.docker.com/engine/install/) Docker desktop
and login using `sudo docker login`. If utilizing GATK-SV Docker images 
from a private container registry or intending to publish the resulting 
images to a registry, ensure that you are logged in with credentials 
that grant you access to the registry.

<Tabs
 groupId="cr"
 defaultValue="gcr"
 values={[
  { label: 'ACR', value: 'acr', },
  { label: 'GCR', value: 'gcr', }
 ]
}>
 <TabItem value="acr">

 You may follow 
 [this documentation](https://learn.microsoft.com/en-us/azure/container-registry/container-registry-authentication?tabs=azure-cli)
 on setting up Docker authentication to an Azure container registry. 
 </TabItem>
 <TabItem value="gcr">

 You may follow
 [this documentation](https://cloud.google.com/artifact-registry/docs/docker/authentication)
 on setting up Docker authentication to a Google container registry. 

 </TabItem>
</Tabs>

### Checkout codebase {#checkout}

Make sure you are on the `git` branch with the code you want to add 
to the GATK-SV Docker images you are building.

```shell
git fetch origin <branch_name>
git checkout origin/<branch_name>
```

## Build and Publish Docker Images {#build}

All the GATK-SV Dockerfiles are hosted under the directory 
[`gatk-sv/dockerfiles/`](https://github.com/broadinstitute/gatk-sv/tree/main/dockerfiles). 
While you can build the GATK-SV Docker images by following the standard 
[Docker image build procedures](https://docs.docker.com/engine/reference/commandline/image_build/),
that can be challenging due to the nested hierarchy of GATK-SV Docker images.
To simplify the process, we have developed a utility script that streamlines the 
Docker image build process 
([`scripts/docker/build_docker.py`](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py)).

In the following, we will explain how to use the utility script for a simple use-case. 
For more advanced and additional functionalities, please refer to the script's documentation,
which you may access it as the following.

```shell
python scripts/docker/build_docker.py --help 
```


In its basic setup, you can use the following command to **build and publish** a GATK-SV Docker image.

```shell
python scripts/docker/build_docker.py \
    --targets <images> \
    --image-tag <tag> \
    --docker-repo <container registry>
```

The arguments used are explained in the following. 

### Determine which images need to be rebuilt {#targets}

You may follow either of the following practices to determine which images to rebuild.

- **Automatic:**
  The script can automatically determine which Docker images need a rebuild 
  based on a list of changed files and cross-referencing them with the 
  table in [this section](/docs/advanced/docker/images#list). 
  Specifically, it takes two git commit SHAs as input, uses `git diff` 
  to extract the list of changed files, and then cross-referencing them 
  with [this table](/docs/advanced/docker/images#list) to identify the Docker 
  images requiring rebuilding. Details can be found on [this page](/docs/advanced/docker/deploy/incremental.md).
  To use this feature, commit the changes first, identify `BASE_SHA` and `HEAD_SHA` using `git log` or GitHub 
  (details on [this page](/docs/advanced/docker/deploy/incremental.md)), 
  and then call the script as follows.

  ```shell
  python scripts/docker/build_docker.py \
      --base-git-commit BASE_SHA \
      --current-git-commit HEAD_SHA
  ```

- **Manual: ** 
  You may refer to the table in [this section](/docs/advanced/docker/images#list)
  to determine which Docker images to rebuild based on the changed files.
  For instance, if you modified any of the files under the
  [`gatk-sv/src/svtk/`](https://github.com/broadinstitute/gatk-sv/tree/main/src/svtk)
  directory, you will need to rebuild the `sv-pipeline` Docker image.
  You can set the list of images to rebuild using the `--targets` argument.
  For instance:
   
  ```shell
  python scripts/docker/build_docker.py \
      --targets sv-pipeline
  ```
     
  You may specify multiple images to rebuild by providing a list of their names. 
  For instance, the following command builds the `sv-pipeline` and the `str` Docker images.
   
  ```shell
  python scripts/docker/build_docker.py \
      --targets sv-pipeline str
  ```
     
Please note that `--targets` and `--base-git-commit --current-git-commit` 
options are mutually exclusive. In other words, you can either manually specify 
images to rebuild, or let the script determine them. 
Combining or avoiding both options is not currently supported.

:::info
Following the steps above, the script builds the specified Docker images 
_and all the images derived from them_, ensuring proper propagation of changes through the pipeline. 
If you want to build only the specified images, you would need to add the `--skip-dependent-images` flag.

<details>
<summary>Build a specific or all the affected Docker images?</summary>

<div>
 <ul>
  <li><b>Build a single image</b> in order to develop and test scripts and Docker images.
   When developing a script encapsulated in a GATK-SV Docker image 
   (e.g., scripts under the <a href="https://github.com/broadinstitute/gatk-sv/tree/main/src">src</a> directory), 
   it's essential to re-build the image and test the script in the image before 
   integrating it into the pipeline. Similarly, when testing a specific task 
   within a workflow, it's easier to reproduce the task in the Docker image 
   where the task is executed. Additionally, if you update the scripts, tools, 
   or configuration within a Dockerfile, rebuilding the respective image is 
   necessary for testing purposes. In all these three cases, building a 
   specific Docker image enables better testing and debugging in isolation.
  </li>
  <li><b>Build all the affected images</b>
   in order to ensure the updated scripts and Dockerfiles are properly 
   propagated through the nested hierarchy of GATK-SV Docker images.
   Once you have successfully developed the specific Docker image and its 
   encapsulated scripts, you will need to ensure that the changes are 
   propagated throughout all the GATK-SV Docker images. 
   This is necessary due to the nested hierarchy
   of GATK-SV Docker images, where derived images should be rebuilt 
   whenever any of their base images are updated. This process 
   involves rebuilding all impacted images affected by the 
   changes you made, ensuring the derived images are updated correctly, 
   and integrating them into the pipeline.
  </li>
 </ul>
</div>
</details>
:::


### Image tag {#tag}

[Docker image tags](https://docs.docker.com/engine/reference/commandline/tag/)
are used to distinguish between different builds of the same image. 
You can use any naming convention for your tags. 
GATK-SV docker images use the following template for tags, 
which you may want to adopt, in particular, if you plan to publish 
your images on the GATK-SV container registries.

```
[Date]-[Release Tag]-[Head SHA 8]
```

where `[Date]` is `YYYY-MM-DD` extracted from the time stamp of the last
commit on the feature branch, `[Release Tag]` is extracted from the latest [pre-]release on GitHub, 
and the `[Head SHA 8]` is the first eight letters of the SHA of the 
last commit on the feature branch.

For example:

```
2023-07-28-v0.28.1-beta-e70dfbd7
```

For automatically composing image tags, you may follow the practices 
used in [GATK-SV CI/CD](https://github.com/broadinstitute/gatk-sv/blob/286a87f3bcfc0b8c811ff789776dd0b135f582e9/.github/workflows/sv_pipeline_docker.yml#L85-L109).



### Specify the container registry {#registry}
The built images are stored on your computer. If you are only developing 
or testing locally, there is no need to push them to a container registry. 
In this case you can avoid providing `--docker-repo <registry>`.

You need to push the images to a container registry if you want to:

- Use the updated Docker images for WDL testing or development;
- Store them on a container registry other than those maintained by the GATK-SV team.
   
The script automatically pushes Docker images to a container registry. 
To use this feature, you may follow these steps:

1. Ensure you are logged into Docker with credentials granting 
push access to the container registry. Please refer to the 
[Docker](#docker) section for details.


2. Provide the `--docker-repo <registry>` argument, 
replacing `<registry>` with the name of your container registry. 
For Google Container Registry (GCR) and Azure Container Registry (ACR), 
the format is generally as follows.

 <Tabs
 groupId="cr"
 defaultValue="gcr"
 values={[
  { label: 'ACR', value: 'acr', },
  { label: 'GCR', value: 'gcr', }
 ]
}>
 <TabItem value="acr">

 Template:

 ```shell
 <registry>.azurecr.io/<repository>/<image>:<tag>
 ```
 
 Example:
 ```shell
 python scripts/docker/build_docker.py \
    --targets sv-pipeline
    --tag v1
    --docker-repo myregistry.azurecr.io/gatk-sv
 ```

 which results in creating the following image:

 ```shell
 myregistry.azurecr.io/gatk-sv/sv-pipeline:v1
 ```

 </TabItem>
 <TabItem value="gcr">

 Template:

 ```shell
 <host name>/<repository>/<image>:<tag>
 ```
 
 Example:
 ```shell
 python scripts/docker/build_docker.py \
    --targets sv-pipeline
    --tag v1
    --docker-repo us.gcr.io/my-repository/gatk-sv
 ```

 which results in creating the following image:

 ```shell
 us.gcr.io/my-repository/gatk-sv/sv-pipeline:v1
 ```

 </TabItem>
</Tabs>



## Post-build

- GATK-SV docker images are mainly intended for use in WDLs. 
  Therefore, it's a good practice to test the newly updated 
  images in related WDLs. This ensures that the updated images function 
  as expected within specific workflows.

- If you were using a Linux VM to build the Docker images, 
  ensure you either stop or delete the VM after building the images. 
  Stopping the VM won't delete the disk, and you'll continue to 
  incur disk usage charges. If you don't want to incur disk costs, 
  you can delete the VM along with all its associated resources. 
  Stopping is preferred over deleting if you intend to reuse the VM.
