---
title: Automated Deployment
description: Build and Publish Images
sidebar_position: 2
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

In the GATK-SV pipeline, the Docker images undergo automated 
processes for building, testing, and publishing as part of 
the CI/CD workflow. These automated procedures guarantee 
that all images are consistently and reproducibly built 
within a standardized Linux VM environment 
(specifically, GitHub Actions). 
This ensures uniformity across all GATK-SV Docker images 
and keeps them synchronized with the latest code-base.


The automated CI/CD pipeline also includes continuous 
testing and regression identification during pull requests. 
This proactive approach allows for the detection and 
resolution of any issues related to image changes or content 
before merging the pull request. 
Consequently, it ensures the reliability and consistency 
of the Docker images, simplifies the review process, 
and maintains the high quality of the pipeline.


Additionally, the automated CI/CD workflow ensures that 
the Docker images are correctly mirrored on multiple 
container registries, specifically Azure Container Registry (ACR) 
and Google Cloud Container Registry (GCR). 
This redundancy guarantees availability and accessibility 
of the images across different platforms.


Latest Docker images are listed in the files, 
with detailed automated deployment descriptions in the following sections.

<Tabs
 groupId="cr"
 defaultValue="gcr"
 values={[
  { label: 'ACR', value: 'acr', },
  { label: 'GCR', value: 'gcr', }
 ]
}>
 <TabItem value="acr">

 ```shell
 gatk_sv_codebase/inputs/values/dockers_azure.json
 ```

 </TabItem>
 <TabItem value="gcr">

 ```shell
 gatk_sv_codebase/inputs/values/dockers.json
 ```

 </TabItem>
</Tabs>


## Workflow Layout

The CI/CD workflow for building, testing, and publishing GATK-SV Docker images 
is defined in [`sv_pipeline.yml`](https://github.com/broadinstitute/gatk-sv/blob/main/.github/workflows/sv_pipeline_docker.yml). 
The [`build_docker.py`](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py) 
script is utilized for building and publishing the images. 
When a pull request is issued against the repository, the images are built, 
and upon merging the pull request, they are published to ACR and GCR.



The workflow consists of three 
[_jobs_](https://docs.github.com/en/actions/learn-github-actions/workflow-syntax-for-github-actions#jobs) 
discussed in the following sections.


### Determine Build Args {#args}
This job is responsible for determining the arguments to be used by the 
`build_docker.py` script, specifically:

- **Determining commit SHAs**:
  Considering the size and number of GATK-SV Docker images, 
  the workflow focuses on building and publishing only the 
  Docker images that are affected by the changes introduced 
  in a pull request (PR).
  You may refer to [this page](/docs/docker/deploy/incremental) 
  on details regarding the incremental build strategy.
  This job determines the commit SHAs of `HEAD` and `BASE`
  commits.

- **Compose image tag**:
  GATK-SV Docker images are tagged with a consistent template 
  to simplify referencing and usage in the pipeline. 
  The tag composition step follows the following template.
  
  ```
  [DATE]-[RELEASE_TAG]-[HEAD_SHA_8]
  ```
  where `[DATE]` represents the `YYYY-MM-DD` format extracted 
  from the timestamp of the last commit on the branch associated 
  with the pull request. `RELEASE_TAG` is extracted from the
  latest [pre-]release on GitHub.
  Additionally, `HEAD_SHA_8` denotes the first eight characters 
  of the `HEAD` commit SHA. The following is an example tag generated
  in this step.
  
  ```
  2023-05-24-v0.27.3-beta-1796b665
  ```
  

### Testing Docker Image Build {#build}
 
The `Test Images Build` job is triggered when a commit is pushed to 
the pull request branch. It is responsible for 
building the Docker images identified by the 
[`Determine Build Args`](#args) 
job. If the Docker image building process fails, 
this job will also fail. The Docker images created 
by this job are not published to GCR or ACR and 
are discarded once the job is successfully completed. 
This job primarily serves for testing purposes during 
the review process, ensuring that the affected images 
can be successfully built and that the changes introduced 
in the pull request do not disrupt the Docker build process.


### Publishing Docker Images {#publish}

The `Publish` job is triggered when a pull request 
is merged or a commit is pushed to the `main` branch. 
Similar to the [`Test Images Build`](#build) job, 
it builds Docker images; however, in addition, 
this job also pushes the built images to the GCR and ACR, 
and updates the list of published images. Specifically, 
this job runs the following steps. 
  

- **Login to ACR**: 
  To authorize access to the Azure Container Registry (ACR), 
  this job logs in to Docker by assuming an Azure service principal.
  The credentials required for the login are defined as 
  [encrypted environment secrets](docs.github.com/en/actions/security-guides/encrypted-secrets).

- **Login to GCR**:
  Similar to ACR, to authorize access to GCR, 
  this job assumes a Google Cloud Platform service account. 
  The secrets related to the service account are defined as 
  [encrypted environment secrets](docs.github.com/en/actions/security-guides/encrypted-secrets).

- **Build and publish to ACR and GCR**:
  Similar to the [build job](#build), this job builds Docker images 
  based on the list of changed files specified using the 
  `HEAD` and `BASE` commit SHA. Additionally, it pushes the 
  built images to both ACR and GCR. It's important to note 
  that the job doesn't rebuild separate images for each registry. 
  Instead, it labels a single image for both ACR and GCR, 
  resulting in an identical image with the same tag and Docker 
  image hash being pushed to both registries. 
  This job will fail if the build or push process encounters any issues.

- **Update the list of published images**:
  GATK-SV maintains two JSON files that store the latest Docker 
  images built and pushed to ACR and GCR. 
  These files require updates whenever a new image is successfully 
  built and published. The `build_docker` script handles the 
  update of the JSON files by adding the latest built and 
  published Docker images for ACR and GCR.

  However, it's important to note that the updated JSON 
  files reside in the GitHub Actions virtual machines, 
  and they are discarded once the GitHub Actions job is 
  completed successfully. To preserve these changes, 
  we need to commit them to the `main` branch from within the 
  GitHub Actions VM as part of the CI/CD process.
  To achieve this, we utilize a dedicated _bot_ account. 
  The steps necessary to perform this are explained  
  in the following.

  - **Login to git using the bot's Personal Access Token (PAT)**:
    This step is necessary to enable the _bot_ account to 
    commit the modified JSON files to the `main` branch 
    and to authorize the _bot_ to push the changes from 
    the GitHub Actions VM to the `main` branch using its credentials.
  
  - **Commit changes and push to the `main` branch**:
    This step configures the Git installation in the 
    GitHub Actions VMs using the _bot_'s credentials. 
    It commits the modified JSON files, which contain 
    the latest built and pushed images. The commit message 
    references the Git commit that triggered the [publish](#publish) job, 
    providing improved tracking of changes in the Git history. 
    Finally, it pushes the commit to the main branch. 
    It's worth noting that Git is intelligent enough 
    to recognize that this push is made from a GitHub 
    Actions environment, preventing it from triggering 
    another publish job. This avoids the issue of 
    infinite triggers of the publish job.

