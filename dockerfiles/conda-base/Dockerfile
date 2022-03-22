ARG CONDA_RELEASE="4.6.14"
ARG PYTHON_RELEASE="3.6.5"
ARG UBUNTU_RELEASE="18.04"
ARG PYTHON_ENV_PATH="/gatk-sv-python"

FROM ubuntu:$UBUNTU_RELEASE as build_conda

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
ARG CONDA_DEP="software-properties-common zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev g++ gfortran ${CONDA_DEP_TRANSIENT}"
# versions of bedtools > 2.27.0 seem to have lost the ability to read gzipped files
# pandas 1.0.0 causes problem with bedtools in aggregate.py
ARG PYTHON_PKGS="pip=21.2.2 setuptools=52.0.0 wheel=0.34.2 bzip2=1.0.8 cython=0.29.14 numpy=1.18.1 pandas=0.25.3 scikit-learn=0.22.1 scipy=1.4.1 intervaltree=3.0.2 matplotlib=3.1.3 natsort=7.0.1 bedtools=2.27.0 pybedtools=0.8.1 pysam=0.14.1=py36_htslib1.7_0 conda-pack=0.7.0"
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR
ARG CONDA_RELEASE
ARG PYTHON_RELEASE
ARG PYTHON_ENV_PATH
ARG CONDA_INSTALL_DIR=/opt/conda
ARG CONDA_BIN=$CONDA_INSTALL_DIR/bin
ARG DEBIAN_FRONTEND=noninteractive
ARG CONDA_PACK_FILE=/tmp/env.tar

ENV PATH=$CONDA_BIN:$PATH

# install conda, create conda environment
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${CONDA_DEP} && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR} && \
    conda install -qy \
          -c conda-forge \
          -c bioconda \
          python=$PYTHON_RELEASE \
          ${PYTHON_PKGS}

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

ENV PATH=$PYTHON_ENV_PATH/bin:$CONDA_BIN:$PATH

# Use conda-pack to move conda packages to a new standalone environment
RUN conda-pack --n-threads -1 -o $CONDA_PACK_FILE && \
    mkdir $PYTHON_ENV_PATH && \
    cd $PYTHON_ENV_PATH && \
    tar xf $CONDA_PACK_FILE && \
    rm $CONDA_PACK_FILE && \
    conda-unpack

# Copy python packages from build_conda. Effectively removes everything but the minimal environment and flattens docker
# layers in the final image.
# NOTE: if using conda-pack, you can no longer use conda (although pip should work fine). This virtual environment is
# "active" automatically as long as $PYTHON_ENV_PATH/bin is in $PATH
FROM ubuntu:$UBUNTU_RELEASE
ARG PYTHON_ENV_PATH
COPY --from=build_conda $PYTHON_ENV_PATH $PYTHON_ENV_PATH
ENV PATH=$PYTHON_ENV_PATH/bin:$PATH
