---
title: Manual Deployment
description: Build and Publish Images
sidebar_position: 3
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';


If you contribute to the GATK-SV codebase, we recommend you build Docker images locally 
to ensure they build successfully and function as intended. The process involves two steps:

- **Build**: Create Docker images from Dockerfiles and store them on your computer.

- **Publish**: Upload the built Docker images to container registries 
(e.g., Google Container registry, or Azure container registry) 
to make them available for use in Terra or Cromwell.
_You may skip this step unless you would like to host the images you built on your own container registry._

To streamline the process, we have developed a 
[script](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py)
that automates both the build and publish steps. 
This section provides guidelines on setting up the environment and running the 
script with a minimal example. For a complete description of the script 
and its various arguments, refer to its help page, as shown in the following.

```shell
python scripts/docker/build_docker.py --help
```

:::danger Linux Machine Required

Only Linux machines (dedicated or virtual) are supported for building GATK-SV Docker images. 
Images created on non-Linux machines (e.g., Apple M1) may not work with Terra or Cromwell execution environment.
The instructions provided on this page assume you are using a Linux Ubuntu machine.
:::



## Setup

### Runtime environment {#runtime}

You may follow the steps in the 
[GCP](https://cloud.google.com/compute/docs/instances/create-start-instance#publicimage)
or [Azure](https://learn.microsoft.com/en-us/azure/virtual-machines/windows/quick-create-portal)
documentation to create a virtual machine (VM) on Google Cloud Platform (GCP) or Microsoft Azure respectively. 
Make sure the VM is built using an Ubuntu image, has at least 8 GB RAM, and some additional 
disk space (e.g., 50 GB should be sufficient).


### Docker {#docker}

[Install](https://docs.docker.com/engine/install/) Docker desktop
and login using `sudo docker login`. 
If you are pulling images from a private container registry,
or intending to publish the resulting images to a registry, 
make sure you login with credentials that grants you with sufficient authorization.

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

```shell
git fetch origin <branch_name>
git checkout origin/<branch_name>
```

## Build and Publish Docker Images {#build}

In its minimal setup, you may use the following command to **build and publish** GATK-SV Docker images.

```shell
python scripts/docker/build_docker.py \
    --targets <images> \
    --image-tag <tag> \
    --docker-repo <container registry>
```

The arguments are explained in the following.


### Determine which images need to be rebuilt {#targets}

You may follow either of the following practices to determine which images to rebuild.

- **Automatic:**
  You may refer to [this page](./images#incremental) for details on this method.
  Briefly, commit the changes first, identify `BASE_SHA` and `HEAD_SHA` using `git log` or GitHub
  and then call the script as follows.

  ```shell
  python scripts/docker/build_docker.py \
      --base-git-commit BASE_SHA \
      --current-git-commit HEAD_SHA
  ```

- **Manual:** 
  You may refer to the table in [this section](./dependencies#list)
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
images to rebuild, or let the script determine them automatically using commit SHAs. 
Combining or avoiding both options is not currently supported.

:::info
Following the steps above, the script builds the specified Docker images 
_and all the images derived from them_. 
You may add the `--skip-dependent-images` flag to build only the explicitly specified images.
:::


### Image tag {#tag}
 
You can use any naming convention for your [tags](https://docs.docker.com/engine/reference/commandline/tag/). 
GATK-SV Docker images are tagged using the following template 
(you may refer to [this section](./automated#args) for details).

```
[Date]-[Release Tag]-[Head SHA 8]
```

For example:

```
2023-07-28-v0.28.1-beta-e70dfbd7
```


### Specify the container registry {#registry}

You may skip this section if you are only developing 
or testing locally; in this case you can avoid providing `--docker-repo <registry>`.
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

Please note that we are currently using GCR, but it has been migrated to Google Artifact Registry.



## Post-build

- GATK-SV docker images are mainly intended for use in WDLs. 
  Therefore, it's a good practice to run the related WDLs with 
  updated images to assert if the images function as expected.

- If you were using a Linux VM to build the Docker images, 
  ensure you either stop or delete the VM after building the images. 
  Stopping the VM won't delete the disk, and you may continue to 
  incur disk usage charges. If you plan on re-using the VM,
  stopping is prefered as it preserves the configuration; 
  otherwise, you may delete the VM and all the associated resources 
  (attached disks in particular).
