# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, sv-pipeline-qc, and sv-pipeline-rdtest
ARG SVBASE_IMAGE=sv-base:latest
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-base-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

FROM $SVBASE_IMAGE

# ##################################################
# # Copy updated virtual environment from VIRTUAL_ENV_IMAGE
ARG VIRTUAL_ENV_PATH
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH

# Add packages needed to support the functioning of the conda environment
ARG CONDA_DEP="software-properties-common zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev g++ gfortran"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${CONDA_DEP} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

##################################################
# copy needed resources from git repo
COPY src/sv-pipeline /opt/sv-pipeline
COPY src/WGD /opt/WGD
COPY src/RdTest /opt/RdTest
COPY src/svtk /opt/svtk
COPY src/svtest /opt/svtest
COPY src/svqc /opt/svqc
ENV PATH="/opt/WGD/bin:${PATH}"
# Install svtk, svtest, svqc and test that they run
# IMPORTANT NOTE: the conda virtual environment provides a special version of pysam
RUN python --version && echo "$PATH" &&\
    echo "about to install svtk" && \
    cd /opt/svtk && \
    pip install -e . && \
    svtk -h && \
    echo "about to install svtest" && \
    cd /opt/svtest && \
    pip install -e . && \
    svtest -h && \
    echo "about to install svqc" && \
    cd /opt/svqc && \
    pip install -e . && \
    svqc -h
