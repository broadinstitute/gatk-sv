---
title: Manual Deployment
description: Build and Publish Images
sidebar_position: 3
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';


If you contribute to the GATK-SV codebase, we recommend you ensure that affected Docker images build successfully and function as intended. The process involves two steps:

1. **Build**: Create Docker images from Dockerfiles.

2. **Publish**: Upload the built Docker images to container registries 
(e.g., Google or Azure container registries, GCR and ACR, respectively) 
to make them available for use in Terra or Cromwell.
_You may skip this step unless you would like to host the images you built on your own container registry._

To streamline the process, we have developed a 
[script](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py)
that automates both the build and publish steps. 
This section provides guidelines on setting up the environment and running the 
script with a minimal example.


:::danger x86 Linux Machine Required
Only Linux machines (dedicated or virtual) are supported for building GATK-SV Docker images. 
In addition, images created on non-Intel processor architectures (e.g., Apple M1) may not function as intended, 
even if the build process runs successfully. 
:::


## Setup an Ubuntu VM

This section outlines steps to follow in order to 
create and connect to a Linux virtual machine (VM)
on a cloud service provider.
You may [skip to the next section](#checkout) if you are using a dedicated Linux machine 
(e.g., a laptop running Ubuntu).


#### 1. Set environment variables

<Tabs
 groupId="cloud"
 defaultValue="gcp"
 values={[
  { label: 'GCP', value: 'gcp', }
 ]
}>
 <TabItem value="gcp">
    ```bash
    export PROJECT_ID="<GOOGLE PROJECT ID>"
    export ZONE_ID="<ZONE ID>"
    
    # Make sure no machine with the following name exist,
    # and you follow VM naming conventions, e.g., all lower-case characters.
    export INSTANCE_NAMES="<VM NAME>"
    ```

 </TabItem>
</Tabs>


#### 2. Create an Ubuntu VM 
You may [skip to the next step](#connect-to-vm) if you have already created a VM.

<Tabs
 groupId="cloud"
 defaultValue="gcp"
 values={[
  { label: 'GCP', value: 'gcp', }
 ]
}>
 <TabItem value="gcp">
    ```bash
    gcloud compute instances create $INSTANCE_NAMES \
      --project=$PROJECT_ID \
      --zone=$ZONE_ID \
      --machine-type=e2-standard-2 \
      --create-disk=auto-delete=yes,boot=yes,device-name=$INSTANCE_NAMES,image=projects/ubuntu-os-cloud/global/images/ubuntu-2310-mantic-amd64-v20240213,mode=rw,size=100
    ```
    Note that this command creates a VM with `100 GiB` disk size, 
    to accommodate for the disk space requirements of GATK-SV Docker images.  

    You may follow the documentation on 
    [this page](https://cloud.google.com/compute/docs/instances/create-start-instance#publicimage)
    for more details on creating a virtual machine on GCP.
 </TabItem>
</Tabs>

:::tip 
The firewall rules of your institute may require you to be on-site or connected 
to the institute's VPN before you can access the cloud resources billed to your institute.
:::


#### 3. Connect to the VM {#connect-to-vm}


<Tabs
 groupId="cloud"
 defaultValue="gcp"
 values={[
  { label: 'GCP', value: 'gcp', }
 ]
}>
 <TabItem value="gcp">
    ```bash
    gcloud compute ssh $INSTANCE_NAMES --project $PROJECT_ID
    ```
    Follow the on-screen prompts for authorizing access to `ssh` credentials.

    <details>
    <summary>Errors running this command</summary>
    <div>
        If you are getting any of the following error messages when you try 
        to connect to the VM immediately after you have created it, 
        it may indicate that the VM is not ready yet, and you may need to
        wait a few minutes before retrying.

        ```bash
        ssh: connect to host [IP address] port 22: Connection refused
        ```

        ```bash
        ERROR: (gcloud.compute.ssh) [/usr/bin/ssh] exited with return code [255].
        username@[IP address]: Permission denied (publickey).
        ```
    </div>
    </details>
 </TabItem>
</Tabs>

#### 4. Install Docker {#docker}
You may [skip to the next step](#checkout) if you have already installed and configured Docker on this VM.

1. Install pre-requisites
    ```bash
    sudo apt-get update && \
    sudo apt-get install ca-certificates curl && \
    sudo install -m 0755 -d /etc/apt/keyrings && \
    sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc && \
    sudo chmod a+r /etc/apt/keyrings/docker.asc && \
    echo \
      "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
      $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
      sudo tee /etc/apt/sources.list.d/docker.list > /dev/null && \
    sudo apt-get update
    ```

2. Install Docker

    ```bash
    sudo apt-get -y install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin && \
    sudo usermod -aG docker ${USER} && \
    newgrp docker
    ```

    You may follow [Docker documentation](https://docs.docker.com/engine/install/ubuntu/) 
    on details on installed Docker on Ubuntu.


3. Login to Docker

    <Tabs
     groupId="cloud"
     defaultValue="gcp"
     values={[
      { label: 'GCP', value: 'gcp', }
     ]
    }>
     <TabItem value="gcp">
        - Run the following command on the VM.
          ```bash
          gcloud auth login
          ```

        - Follow the on-screen prompts, it will display a URL that you need to copy-paste it
          on the browser of your computer (_not_ the VM).

        - Follow the prompts on your browser, and login with an account that will provide you
          with access to the GCR repository. If you are planning on _publishing_ images you 
          build to GCR, you need to make sure you account has [sufficient access](https://cloud.google.com/artifact-registry/docs/docker/pushing-and-pulling#required_roles)
          to GCR. 

        - Configure Docker with your credentials.

          ```bash
          gcloud auth configure-docker
          ```

        You may refer to [this page](https://cloud.google.com/artifact-registry/docs/docker/authentication)
        for more details on configure Docker to access GCR.
     </TabItem>
    </Tabs>


## Checkout codebase {#checkout}

1. Clone the repository or its fork that contains the branch with the changes 
that you want to build the Docker images based-off. 

    ```shell
    git clone https://github.com/broadinstitute/gatk-sv && cd gatk-sv
    ```

2. Checkout the branch containing your changes.

    ```shell
    git checkout <BRANCH_NAME>
    ```

## Build and Publish Docker Images {#build}

In its minimal setup, you may use the following command to **build and publish** GATK-SV Docker images.

```shell
python3 scripts/docker/build_docker.py \
    --targets <IMAGES> \
    --image-tag <TAG> \
    --docker-repo <CONTAINER_REGISTRY>
```

The arguments are explained in the following.

- [`--targets`](#targets)
- [`--image-tag`](#tag)
- [`--docker-repo`](#registry)

### `--targets` {#targets}

You may follow either of the following approaches to determine which images to rebuild.

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

- **Automatic (advanced):**
  You may refer to [this page](./images#incremental) for details on this method.
  Briefly, you may take the following steps. 

  1. `git commit` the changes. 
  2. Identify `BASE_SHA` and `HEAD_SHA` using `git log` or GitHub. 
     You may use the following commands to get these SHAs.

     ```shell
     export \
       HEAD_SHA=$(git log -1 --pretty=format:"%H") \
       BASE_SHA=$(git merge-base main $(git branch --show-current))
     ```
     Note that, you may need to [modify these commands](https://git-scm.com/docs/git-merge-base) if your branch has a complicated git history.

  3. Run the script using `--base-git-commit` and `--current-git-commit` instead of `--targets`. 
  ```shell
  python scripts/docker/build_docker.py \
      --base-git-commit <BASE_SHA> \
      --current-git-commit <HEAD_SHA>
  ```
     
Please note that `--targets` and `--base-git-commit --current-git-commit` 
options are mutually exclusive. In other words, you can either manually specify 
images to rebuild, or let the script determine them automatically using commit SHAs; 
combining or avoiding both options is not currently supported.

:::info
Following the steps above, the script builds the specified Docker images 
_and all the images derived from them_. 
You may add the `--skip-dependent-images` flag to build only the explicitly specified images.
:::


### `--image-tag` {#tag}
 
You may use any naming convention for the Docker image 
[tags](https://docs.docker.com/engine/reference/commandline/tag/). 
GATK-SV Docker images are tagged using the following template 
(you may refer to [this section](./automated#args) for details).

```
[Date]-[Release Tag]-[Head SHA 8]
```

For example:

```
--image-tag 2023-07-28-v0.28.1-beta-e70dfbd7
```


### ` --docker-repo` {#registry}

If you are only testing GATK-SV Docker image build, 
you may skip this section and avoid providing `--docker-repo <registry>`.
However, if you need to push image to container registries, 
need images for WDL testing, or need to host the images on a container registry 
other than those maintained by the GATK-SV team. 
   
The `build_docker.py` script automatically pushes Docker images to a container registry
when `--docker-repo <registry>` is provided, replacing `<registry>` with the container registry you want to use.
When providing this argument, ensure that you are logged into Docker with 
credentials granting push access to the registry,
You may configure and set the registry as the following.

 <Tabs
 groupId="cr"
 defaultValue="gcr"
 values={[
  { label: 'ACR', value: 'acr', },
  { label: 'GCR', value: 'gcr', }
 ]
}>
 <TabItem value="acr">

 - You may follow [these steps](https://learn.microsoft.com/en-us/azure/container-registry/container-registry-get-started-portal?tabs=azure-cli)
   if you have not configured a container registry.
 - Once configured, you may set `<registry>` in the following template.

    ```shell
    <REGISTRY>.azurecr.io/<REPOSITORY>/<IMAGE>
    ```
   
    Example:

    ```shell
    myregistry.azurecr.io/gatk-sv
    ```

 </TabItem>
 <TabItem value="gcr">

 - You may follow [these steps](https://cloud.google.com/artifact-registry/docs/repositories/create-repos)
   if you have not configured a container registry.
 - Once configured, you may set `<registry>` in the following template. 

     ```shell
     <HOST_NAME>/<REPOSITORY>/<IMAGE>
     ```
    
    Example:
    ```shell
    us.gcr.io/my-repository/gatk-sv
    ```

 </TabItem>
</Tabs>


## Post-build

- GATK-SV docker images are mainly intended for use in WDLs. 
  Therefore, it's a good practice to run the related WDLs with 
  updated images to assert if the images function as expected.

- If you were using a Linux VM to build the Docker images, 
  ensure you either stop or delete the VM after building the images. 
  Stopping the VM won't delete the disk, and you may continue to 
  incur disk usage charges. If you plan on re-using the VM,
  stopping is preferred as it preserves the configuration; 
  otherwise, you may delete the VM and all the associated resources 
  (attached disks in particular).
