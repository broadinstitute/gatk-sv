ARG VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG CONDA_RELEASE="4.6.14"
ARG PYTHON_RELEASE="3.6.5"

# build changes here, then poof the completed changes into the final image
FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

##################################################
# shared R packages by sv-pipeline* dockers
ARG SV_PIPE_SHARED_R_PKGS="MASS e1071"
ARG VIRTUAL_ENV_PATH
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake \
                 file \
                 g++ \
                 gfortran \
                 liblapack-dev \
                 libopenblas-dev \
                 libxml2-dev && \
    install_deprecated_R_package.sh https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz

##################################################
# conda and all packages
# Steps:
#  1. OS libs
#  2. install miniconda, enable for all users (instead of current user only), update and cleanup
#  3. install some commonly used conda packages
#  4. final clean up
#  5. special note:  Pysam 0.14.1 - Installing from source prevents the following bug:
#                    python: vcf.c:3482: bcf_update_format: Assertion `!fmt->p_free' failed.
#                    /cromwell_root/script: line 31: 22 Aborted (core dumped) /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py
# NB: The conda installation here overrides samtools 1.9 and bcftools 1.9 with versions 1.7 due to
# a pysam dependency.
ARG CONDA_DEP_TRANSIENT="make git wget curl apt-transport-https ca-certificates gnupg-agent software-properties-common"
ARG CONDA_DEP="zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev g++ gfortran ${CONDA_DEP_TRANSIENT}"
# versions of bedtools > 2.27.0 seem to have lost the ability to read gzipped files
# pandas 1.0.0 causes problem with bedtools in aggregate.py
ARG PYTHON_PKGS="pip=21.2.2 setuptools=52.0.0 wheel=0.34.2 bzip2=1.0.8 cython=0.29.14 numpy=1.18.1 pandas=0.25.3 scikit-learn=0.22.1 scipy=1.4.1 intervaltree=3.0.2 matplotlib=3.1.3 natsort=7.0.1 bedtools=2.27.0 pybedtools=0.8.1 pysam=0.14.1=py36_htslib1.7_0 conda-pack=0.7.0"
ARG SV_BASE_R_PKGS="optparse BiocManager"
ARG CONDA_RELEASE
ARG PYTHON_RELEASE
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR
ARG CONDA_RELEASE
ARG DEBIAN_FRONTEND=noninteractive

# get prerequisites for conda
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${CONDA_DEP}
# install conda
ARG CONDA_INSTALL_DIR=/opt/conda
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR}

# create conda environment
ARG CONDA_BIN=$CONDA_INSTALL_DIR/bin
ENV PATH=$CONDA_BIN:$PATH
ARG PYTHON_RELEASE
RUN conda install -qy -c conda-forge python=$PYTHON_RELEASE mamba

# install python packages
RUN mamba install -qy \
          -c conda-forge \
          -c bioconda \
          conda-pack=0.7.0 \
          $PYTHON_PKGS

# remove pysam as installed by conda, and replace it with special version to fix bug in svtk
RUN conda remove --force pysam && \
    git clone https://github.com/pysam-developers/pysam.git /opt/pysam && \
    cd /opt/pysam && \
    git checkout af6ff23322af451b8ab6ea78b53676c45ae0a779 && \
    pip install --no-cache-dir .

# NOTE: if you were going to copy over the conda environment, you'd want to clean it out like so:
#RUN conda clean -ay --force-pkgs-dirs
# However, conda-pack cleans up the environment itself, and is actually confused by running conda clean first. So this
# comment block is left as a warning to NOT do this.

ENV PATH=$CONDA_BIN:$PATH

# Use conda-pack to move conda packages to a new standalone environment
ARG CONDA_PACK_FILE=/tmp/env.tar
RUN conda-pack --n-threads -1 -o $CONDA_PACK_FILE && \
    mkdir -p $VIRTUAL_ENV_PATH && \
    cd $VIRTUAL_ENV_PATH && \
    tar xf $CONDA_PACK_FILE && \
    rm $CONDA_PACK_FILE && \
    conda-unpack

#########################################################################################
# Copy packages from build_virtual_env. Poofs the changes to the virtual environment without any extra cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
