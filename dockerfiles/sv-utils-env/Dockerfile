ARG UBUNTU_RELEASE="22.04"
ARG PYTHON_VIRTUAL_ENV_IMAGE=samtools-cloud-virtual-env:latest
ARG CONDA_INSTALL_DIR=/opt/conda
ARG CONDA_ENV_NAME="gatk-sv"

############# get python environment and add to it
FROM $PYTHON_VIRTUAL_ENV_IMAGE as python_virtual_env

# install any build dependencies, plus ghostscript (needed for merging PDFs)
ARG BUILD_DEPS="g++ make apt-transport-https ca-certificates wget zlib1g-dev"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends $BUILD_DEPS

# install python packages
RUN mamba install -qy --freeze-installed -n $CONDA_ENV_NAME -c conda-forge -c bioconda \
          psutil=5.8.0 pympler=0.9 tqdm=4.62.3 \
          numpy=1.23.0 pandas=1.4.3 scipy=1.7.3 dill=0.3.5.1 \
          seaborn=0.11.2 matplotlib-base=3.5.1 ghostscript=9.54.0 \
          cython=0.29.28 \
          pytest=7.1.2 pytest-cov=3.0.0

# As of March 29, 2022 the newest pysam version is 0.19.0, however that version does not allow setting END < POS, which
# causes a problem for our BND handling. The newest version that is compatible with our BND processing is 0.15.4
# In general, pysam needs to be pip installed or built, because you want to use the installed HTSLIB rather than copy in
# duplicate libraries. The older pysam that needs to be installed for BND issues has problems with cython during pip
# install, so must be built from source. One last pysam wrinkle: it's install can't handle the newest versions of
# setuptools, so we're going to downgrade then re-upgrade them just for pysam
ARG PYSAM_VERSION=0.15.4
RUN export SETUPTOOLS_VERSION=$(python -c 'import setuptools; print(setuptools.__version__)') && \
    pip install setuptools==57.5.0 && \
    wget -q https://github.com/pysam-developers/pysam/archive/refs/tags/v$PYSAM_VERSION.tar.gz && \
    tar -xzf v$PYSAM_VERSION.tar.gz && \
    cd pysam-$PYSAM_VERSION && \
    python setup.py build && \
    python setup.py install && \
    pip install setuptools==$SETUPTOOLS_VERSION

# clean unneeded stuff
RUN conda clean -ay --force-pkgs-dirs

#################### build changes here, then copy the completed changes into the final image
FROM ubuntu:$UBUNTU_RELEASE
ARG CONDA_INSTALL_DIR
COPY --from=python_virtual_env $CONDA_INSTALL_DIR $CONDA_INSTALL_DIR
ARG CONDA_BIN=$CONDA_INSTALL_DIR/bin
ARG CONDA_ENV_NAME
ENV CONDA_ENV_NAME=$CONDA_ENV_NAME
ARG CONDA_ENV_BIN=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME/bin
ENV PATH=$CONDA_ENV_BIN:$CONDA_BIN:$PATH
