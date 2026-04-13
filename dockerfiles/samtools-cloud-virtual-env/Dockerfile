ARG UBUNTU_RELEASE="22.04"
ARG CONDA_RELEASE="py39_25.3.1-1"  # releases listed here: https://repo.anaconda.com/miniconda/
ARG PYTHON_RELEASE="3.10.4"  # get available releases from conda search python
ARG CONDA_INSTALL_DIR=/opt/conda
ARG CONDA_ENV_NAME="gatk-sv"

#################### build changes here, then copy the completed changes into the final image
FROM ubuntu:$UBUNTU_RELEASE as build

# get prerequisites for installing conda
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 apt-transport-https ca-certificates wget

# install conda
ARG CONDA_RELEASE
ARG CONDA_INSTALL_DIR
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR}

# update python to desired version, and install mamba for faster package management
ARG PYTHON_RELEASE
ARG CONDA_BIN=$CONDA_INSTALL_DIR/bin
ARG CONDA_ENV_BIN=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME/bin
ENV PATH=$CONDA_ENV_BIN:$CONDA_BIN:$PATH

# # create conda environment
ARG CONDA_ENV_NAME
# install the desired python version, and just the packages needed by google-cloud-cli
RUN conda install -qy -c conda-forge mamba=2.1.1 && \
    mamba create -n $CONDA_ENV_NAME -qy -c conda-forge python=$PYTHON_RELEASE \
      pip=25.1.1 setuptools=80.1.0 wheel=0.45.1 crcmod=1.7

# at this stage, we only need python for google-cloud-cli, however, we will want to use conda to install more packages
# later, so copy the whole conda install over, not just the environment
# clean unneeded stuff first, to decrease size
RUN conda clean -ay --force-pkgs-dirs

#########################################################################################
# Copy packages from build. Get the final image without intermediate changes or extra files
FROM ubuntu:$UBUNTU_RELEASE
ARG CONDA_INSTALL_DIR
ENV CONDA_INSTALL_DIR=$CONDA_INSTALL_DIR
COPY --from=build $CONDA_INSTALL_DIR $CONDA_INSTALL_DIR
ARG CONDA_BIN=$CONDA_INSTALL_DIR/bin
ARG CONDA_ENV_NAME
ENV CONDA_ENV_NAME=$CONDA_ENV_NAME
ARG CONDA_ENV_BIN=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME/bin
ENV PATH=$CONDA_ENV_BIN:$CONDA_BIN:$PATH
