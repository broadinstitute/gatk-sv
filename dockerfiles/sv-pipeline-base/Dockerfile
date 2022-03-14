# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, sv-pipeline-qc, and sv-pipeline-rdtest
#
ARG SVBASE_IMAGE=sv-base:latest
ARG CONDABASE_IMAGE=conda-base:latest
ARG PYTHON_ENV_PATH="/gatk-sv-python"

# establish alias to copy from later:
FROM $CONDABASE_IMAGE as conda_base

FROM $SVBASE_IMAGE

##################################################
# shared R packages by sv-pipeline* dockers
ARG SV_PIPE_SHARED_R_PKGS="MASS e1071"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake \
                 file \
                 g++ \
                 gfortran \
                 liblapack-dev \
                 libopenblas-dev \
                 libxml2-dev && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
    /opt/install_deprecated_R_package.sh https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy remove make cmake automake && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# ##################################################
# # Copy packed conda virtual environment from CONDABASE_IMAGE
ARG PYTHON_ENV_PATH
COPY --from=conda_base ${PYTHON_ENV_PATH} ${PYTHON_ENV_PATH}
# Add virtual environment executables to PATH. Now environment is installed and "activated"
ENV PATH="${PYTHON_ENV_PATH}/bin:${PATH}"

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
