# PanGenie docker file.
# In this Dockerfile PanGenie is installed via conda (recommended approach
# by its maintainers https://bitbucket.org/jana_ebler/pangenie/src/master/).
#
# A note on the base image, miniconda3:
# Installation via conda requires activation of PanGenie conda environmet;
# however, that requires initialization bash to use conda (error message:
# "Your shell has not been properly configured to use 'conda activate'".).
# Conda initialization and bash profile reload via the following does not
# resolve the issue.
#   RUN bash miniconda.sh -b -p /opt/conda && \
#       /opt/conda/bin/conda init bash && \
#       bash -c "source ${HOME}/.bashrc"
#
# A similar issue is detailed in the following blog post; therefore, as
# suggested in this blog post, this docker files is using
# 'continuumio/miniconda3' as base image.
#   https://pythonspeed.com/articles/activate-conda-dockerfile/

FROM continuumio/miniconda3:4.10.3p0

RUN apt-get update && apt-get install --no-install-recommends -qqy \
    build-essential \
    git \
    make \
    cmake \
    wget \
    pkg-config

RUN git clone https://jana_ebler@bitbucket.org/jana_ebler/pangenie.git && \
    cd pangenie && \
    # Since no releases on bitbucket, pinning this to it current latest
    # commit on the master branch.
    git checkout 7c9a8b7eb9c5ff03cbf8a7ef128e86241d259a76 && \
    conda env create -f environment.yml

# Make RUN commands use the new environment:
RUN echo "conda activate pangenie" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Install PanGenie
# Without setting `PKG_CONFIG_PATH`, cmake/make cannot resolve
# jellyfish dependency.
# The jellyfish version is pinned to its current latest release.
RUN conda install -y jellyfish=2.2.10 -c bioconda && \
    export PKG_CONFIG_PATH="/opt/conda/envs/pangenie/lib/pkgconfig/" && \
    cd pangenie && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make
ENV PATH="/pangenie/build/src/:$PATH"

RUN apt-get -qqy clean && \
    rm -rf  /tmp/* \
            /var/tmp/* \
            /var/cache/apt/* \
            /var/lib/apt/lists/* \
            /usr/share/man/?? \
            /usr/share/man/??_*
