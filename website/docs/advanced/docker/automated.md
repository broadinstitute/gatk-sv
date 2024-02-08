---
title: Automated Deployment
description: Build and Publish Images
sidebar_position: 2
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';


GATK-SV Docker images are automatically built, tested, and hosted on 
container registries. The automation is part of the pipeline’s continuous 
integration and continuous delivery (CI/CD), which ensures the 
images are built and tested consistently and reproducibly in standardized Linux virtual machines.


The automation pipeline runs on GitHub Actions; it performs a regression 
test as part of every pull request. When a pull request is merged, the automation 
pipeline publishes images on the Google Container Registry (GCR) 
and Azure Container Registry (ACR), and updates their references.


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

 
 [gatk-sv/inputs/values/dockers_azure.json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers_azure.json)
 

 </TabItem>
 <TabItem value="gcr">

 [gatk-sv/inputs/values/dockers.json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json)
 

 </TabItem>
</Tabs>


:::info
The detailed explanation of the automation workflow provided on this page
is intended for users who need to configure the CI/CD workflow on
their own fork of GATK-SV's GitHub repository to host Docker images on 
their own container registries.


If you only need the list of latest Docker images, you may refer to the above-listed files.
:::


## Workflow Layout

The automation workflow is defined in 
[`sv_pipeline.yml`](https://github.com/broadinstitute/gatk-sv/blob/main/.github/workflows/sv_pipeline_docker.yml), 
and it utilizes the 
[`build_docker.py`](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py)
script for building and publishing Docker images.
The workflow consists of three 
[_jobs_](https://docs.github.com/en/actions/learn-github-actions/workflow-syntax-for-github-actions#jobs) 
discussed in the following sections.

- [Determine build arguments](#args);
- [Regression testing](#build) (pull request and merge);
- [Publishing Docker images](#publish) (merge only).

### Determine Build Args {#args}
This job is responsible for determining the arguments to be used by the 
`build_docker.py` script, specifically:

- **Determining commit SHAs**:
  Considering the size and number of GATK-SV Docker images, 
  the workflow builds and publishes only the 
  Docker images affected by the changes introduced 
  in a pull request.
  You may refer to [this page](/docs/advanced/docker/deploy/incremental) 
  on details regarding the incremental build strategy.
  This job determines the commit SHAs of `HEAD` and `BASE`
  commits.

- **Compose image tag**:
  The images are tagged with a consistent template as the following.
  
  ```
  [DATE]-[RELEASE_TAG]-[HEAD_SHA_8]
  ```

  - `[DATE]` is in `YYYY-MM-DD`, and is extracted 
  from the timestamp of the last commit on the branch associated 
  with the pull request. 
  - `RELEASE_TAG` is extracted from the
  latest [pre-]release on GitHub.
  - `HEAD_SHA_8` denotes the first eight characters 
  of the `HEAD` commit SHA. 
  
  The following is an example of a tag generated
  in this step.
  
  ```
  2023-05-24-v0.27.3-beta-1796b665
  ```
  

### Testing Docker Image Build {#build}

This job is triggered when **a commit is pushed to the pull request branch.**
It serves the purpose of regression testing of the Docker images.
It builds Docker images according to the arguments determined in [`Determine Build Args`](#args). 
This job fails if it fails to build the Docker images, 
and the images it builds are discarded when the job finishes.


### Publishing Docker Images {#publish}

This job is triggered when **a pull request is merged or a commit is pushed to the `main` branch.**
Similar to the [`Test Images Build`](#build) job, 
it builds Docker images; however, in addition, 
this job also pushes the built images to the GCR and ACR, 
and updates their list. 
The job fails if it cannot successfully run all the steps.
Specifically, it takes the following steps.
  

- **Login**
  to container registries in order to push the built images.
  The job obtains authorization to push to Google and Azure container registries 
  by assuming a Google service account and an Azure service principal, respectively. 
  The credentials required to assume these identities are defined as 
  [encrypted environment secrets](https://docs.github.com/en/actions/security-guides/encrypted-secrets).


- **Build and publish to ACR and GCR**:
  Similar to the [build job](#build), this job builds Docker images 
  based on the list of changed files specified using the 
  `HEAD` and `BASE` commit SHA. However, in addition, it pushes the 
  built images to both ACR and GCR. It's important to note 
  that the images pushed to GCR and ACR are identical, and they only differ in their tag.

- **Update the list of published images**:
  Once the newly built images are successfully pushed, 
  this job updates the JSON files containing the list of images (i.e., `dockers*.json`),
  and commit and push the changes to the `main` branch.
  To achieve this, we use a _bot_ account; specifically:

  - Login to git using the bot's Personal Access Token (PAT)
    in order to authorize it to push to the `main` branch.
  
  - Configure the Git installation in the GitHub Actions VMs using the _bot_'s credentials. 

  - Commit the changed files. The commit message references the 
    Git commit that triggered the [publish](#publish) job. 
  
  - Push the commit to the main branch.
  
  It is worth noting that GitHub recognizes that this push to the `main` branch is made from a GitHub 
  Actions environment, hence it does not trigger another [Publish](#publish) job,
  avoiding infinite triggers of this job.
