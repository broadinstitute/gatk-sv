---
title: Docker Images
description: Building and hosting of the Docker Images
sidebar_position: 4
---

The tools needed to run the GATK-SV pipeline are organized in 
different Docker images. These images are automatically 
built in GATK-SV [_Continuous Integration and Continuous Delivery_ (CI/CD)](https://github.com/broadinstitute/gatk-sv/tree/main/.github/workflows)
and hosted on Google Container Registry (GCR). 
The [`dockers.json`](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json)
file contains an up-to-date list of Docker images.


:::info
The process of building Docker images, uploading to a container registry, and updating the `docker.json`
file is completely automated. You do not need to re-run this process unless for development purposes.

The GATK-SV CI/CD does not update the references to Docker images in your workflows. 
Therefore, you need to consult with `dockers.json` for the latest current Docker images to be used 
in your workflows. 
:::

### When should a Docker image be updated?

A Docker image should be updated if any of the files included in 
the image is modified. In general, if you need to make a change to a script 
(basically anything in the pipeline that’s not a WDL), you will need to build 
a new docker image and update the parameters in the test JSON files in order 
to propagate your change. 

You may take the following steps to build GATK-SV Docker images.

1. You should have already tested your changes by running just the 
   updated script inside the old Docker image, because it is 
   faster and easier to debug interactively rather than building the image, 
   submitting a job, then scouring the log files for any errors.

2. You have two options to build docker images: run `build_docker.py` use `Docker build`.

   * The GATK-SV docker images have hierarchical dependencies that 
     could be non-trivial to manually resolve. Hence we developed a
     utility script that streamlines building the images. You may run 
     the script as the following.

     ```shell
     $ python scripts/docker/build_docker.py --targets LIST_OF_IMAGES_TO_REBUILD
     ```
   
     Please refer to the help (i.e., `python build_docker.py --help`) or the 
     [README](https://github.com/broadinstitute/gatk-sv/tree/main/scripts/docker)
     for more details and other arguments of this script.

   * To manually build a Docker image you may run the following, which is a general 
     way of building any Docker image. However, due to the hierarchical dependencies
     between GATK-SV Docker images, you would need to make sure all the dependant 
     Docker images are re-built. 

     ```shell
     $ docker image build -f dockerfiles/[image directory]/Dockerfile -–no-cache -–tag [tag] .
     ```
     
     Note the `.` at the end of the command. It is required to invoke `docker image build` 
     in this way since some dockerfiles copy scripts from gatk-sv to the docker image, 
     and they could correctly resolve path if invoked this way. 

4. After you build your docker image, test the WDL you altered (or perhaps several 
   affected modules) by changing the docker image argument in the inputs JSON to 
   your new image (the format is `repo/image:tag` from the arguments you provided 
   to the build script).

   - You only need to update the docker argument for the image that you just 
     built (i.e., if you built the target sv-pipeline-qc, you just need to 
     change the sv_pipeline_qc_docker argument(s), not the 
     sv_base_mini_docker or sv_pipeline_docker arguments). 

   - Then, run the WDL(s) with cromshell or Terra with the new docker image to make sure everything works.
