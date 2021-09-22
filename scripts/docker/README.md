## Overview

The scripts here are to be used for building all, or any combination of dockers that are listed in `dockerfiles` directory
(and have their dependencies specified in build_docker.py).

We support both building from your local git repo, and remote Github repo.

### Dependencies
This set of python script depends on Python 3.

Conda packages `termcolor` and `pprint` are assumed to be available.

### Usage

1. Checkout your working branch:
    ```
    $ cd ${REPO_ROOT}
    $ git checkout my_branch
    ```

2. Run the docker build script:
    ```
    $ cd scripts/docker
    $ python build_docker.py --targets sv-pipeline --image-tag my-branch-a3fd92 --dockerhub-root myrepo 
    ```
    This example will build and push `myrepo/sv-pipeline:my-branch-a3fd92`. All required base images will 
    also be included automatically. Additional options can be viewed by running `python build_docker.py`.

### Notes

* Multiple target images can be specified like so: `--targets sv-pipeline wham manta` or `--targets all`

* To build docker(s) that you desire, we recommend you look at `build_docker_template.ipynb`, which consists of example use cases when building from local and remote git resources

* A full build (e.g. using `--targets all`) may take a few hours

* Due to license restrictions, we cannot provide the binaries required for running the MELT docker. Users must obtain a license and download `MELTv2.0.5_patch.tar.gz` to `dockerfiles/melt` in order to build the MELT docker for module 00a. Alternatively, users may run MELT locally, and upload the VCFs to a GCS bucket to use as input to modules 00b and 00c. 

* When building from local files, we cautiously refuse to build when there are uncommitted changes, and/or un-tracked files, unless you specifically turn off that protection with `--disable-git-protect`.

* When building from Github tag/hash values, we assume you'd provide an **existing** temporary locations where we'll pull the files to. When we are done, the pulled files will be cleaned up, but not the temporary location you provided. When you'd like to pull from Github using SSH instead of https, simply turn on flag `--use_ssh`.

* `test_build_docker.ipynb` is for debugging and testing purposes

* Use `update_json_docker.sh` to update dockers tags across all json parameter files
