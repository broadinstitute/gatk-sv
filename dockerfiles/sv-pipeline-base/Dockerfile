# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, sv-pipeline-qc, and sv-pipeline-rdtest

# NB: The conda installation here overrides samtools 1.9 and bcftools 1.9 with versions 1.7 due to
# a pysam dependency.

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG GATKSV_PIPELINE_V1_RESOURCES_IMAGE=gatksv-pipeline-v1-resources:latest
ARG SVBASE_IMAGE=gatksv/sv-base:latest

####################################################################################################
### stage 1 for moving needed files to a promised location
FROM ${GATKSV_PIPELINE_V1_RESOURCES_IMAGE}
RUN mkdir -p /opt/ && \
    mv ${RESOURCES}WGD /opt/WGD && \
    mv ${RESOURCES}svtk /opt/svtk && \
    mv ${RESOURCES}sv-pipeline /opt/sv-pipeline && \
    mv ${RESOURCES}RdTest /opt/RdTest && \
    mv ${RESOURCES}svtest /opt/svtest && \
    mv ${RESOURCES}svqc /opt/svqc

####################################################################################################
### stage 2
FROM ${SVBASE_IMAGE}

##################################################
# shared R packages by sv-pipeline* dockers
ARG SV_PIPE_SHARED_R_PKGS="MASS e1071"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
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
ARG CONDA_DEP_TRANSIENT="make git wget"
ARG CONDA_DEP="software-properties-common zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libblas-dev liblapack-dev libatlas-base-dev g++ gfortran ${CONDA_DEP_TRANSIENT}"
# versions of bedtools > 2.27.0 seem to have lost the ability to read gzipped files
# pandas 1.0.0 causes problem with bedtools in aggregate.py
ARG PYTHON_PKGS="setuptools=52.0.0 wheel=0.34.2 bzip2=1.0.8 cython=0.29.14 numpy=1.18.1 pandas=0.25.3 scikit-learn=0.22.1 scipy=1.4.1 intervaltree=3.0.2 matplotlib=3.1.3 natsort=7.0.1 bedtools=2.27.0 pybedtools=0.8.1 pysam=0.14.1=py36_htslib1.7_0"
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR="/opt/conda"
ARG CONDA_RELEASE="4.6.14"
ARG CONDA_BIN=${CONDA_INSTALL_DIR}/bin
ARG CONDA_CMD=${CONDA_BIN}/conda
ENV PATH=${CONDA_BIN}:$PATH
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${CONDA_DEP} && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR} && \
    ln -s ${CONDA_INSTALL_DIR}/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo "${CONDA_CMD} activate base" >> ${HOME}/.bashrc && \
    /bin/bash -c "source ${HOME}/.bashrc" && \
    conda install -qy \
          -c conda-forge \
          -c bioconda \
          python=3.6.5 \
          ${PYTHON_PKGS} && \
    git clone https://github.com/pysam-developers/pysam.git /opt/pysam && \
    cd /opt/pysam && \
    git checkout af6ff23322af451b8ab6ea78b53676c45ae0a779 && \
    pip install -e . && \
    rm -rf ${HOME}/.cache/pip && \
    conda clean -tipsy && \
    rm /tmp/miniconda.sh && \
    apt-get -qqy remove --purge ${CONDA_DEP_TRANSIENT} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

##################################################
# copy needed resources from previous stage 
# Add sv-pipeline
COPY --from=0 /opt/sv-pipeline /opt/sv-pipeline

# Add WGD
COPY --from=0 /opt/WGD /opt/WGD
ENV PATH="/opt/WGD/bin:${PATH}"

# Add RdTest
COPY --from=0 /opt/RdTest /opt/RdTest

# Install svtk and test that it runs 
# IMPORTANT NOTE: this is using, because of the conda packages installation above, a special version of pysam
COPY --from=0 /opt/svtk /opt/svtk

RUN echo "about to install svtk" && \
    python --version && cd /opt/svtk && \
    pip install -e . && \
    svtk -h

# Install svtest and test that it runs
COPY --from=0 /opt/svtest /opt/svtest

RUN echo "about to install svtest" && \
    python --version && cd /opt/svtest && \
    pip install -e . && \
    svtest -h

# Install svqc and test that it runs
COPY --from=0 /opt/svqc /opt/svqc

RUN echo "about to install svqc" && \
    python --version && cd /opt/svqc && \
    pip install -e . && \
    svqc -h

##################################################
########### TODO BLOCK: EVEN IF THESE FOLLOWING INSTALLATIONS ARE NEEDED, WE SHOULD REALLY MAKE THEM LOCAL

# # Python 3.6 environment
# RUN conda update -qy -n base conda && \
#     conda install -qy -c bioconda python=3.6 numpy scipy cython natsort boto3 pandas scikit-learn

# # Shuang: I'm taking out this statement because when I searched for `gsutil signurl` in the repo, only WDLs in the subdir "sv-pipline" use that, 
# #         and my understanding is that WDLS there will not be used in production.
# # Python 2.7 environment with pyopenssl for `gsutil signurl`
# RUN conda create -y \
#           -n py27 \
#           python=2.7 pyopenssl && \
#     /bin/bash -c "source activate py27 && pip install -q --no-cache-dir google-compute-engine"
# ENV CLOUDSDK_PYTHON /opt/conda/envs/py27/bin/python

# # Shuang: I'mm also taking this out, since I do not think Google cloud functionality is to be installed through pip, 
# #         and kubernetes is just plain useless for us at this stage.
# # Q: what for?, especially kubernetes
# RUN pip install -q --no-cache-dir \
#                 google-cloud-storage \
#                 kubernetes

# # Shuang: the specific fix related to the sv-pipline issue is moved up, to avoid installing packages twice
# # Pysam 0.14.1 - Installing from source prevents the following bug:
# # python: vcf.c:3482: bcf_update_format: Assertion `!fmt->p_free' failed.
# # /cromwell_root/script: line 31: 22 Aborted (core dumped) /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py
# ARG DEBIAN_FRONTEND=noninteractive
# RUN apt-get -qqy update --fix-missing && \
#     apt-get -qqy dist-upgrade && \
#     apt-get -qqy install --no-install-recommends \
#                  cmake \
#                  git \
#                  make \
#                  libcurl4-openssl-dev && \
#     conda install -qy \
#           -c conda-forge \
#           -c bioconda \
#           bzip2 \
#           cython \
#           pybedtools \
#           pysam=0.14.1=py36_htslib1.7_0 && \
#     git clone https://github.com/pysam-developers/pysam.git /opt/pysam && \
#     cd /opt/pysam && \
#     git checkout af6ff23322af451b8ab6ea78b53676c45ae0a779 && \
#     pip install -e . && \
#     conda clean -tipsy && \
#     apt-get -qqy remove --purge git make cmake && \
#     apt-get -qqy autoremove --purge && \
#     apt-get -qqy clean && \
#     rm -rf /tmp/* \
#            /var/tmp/* \
#            /var/cache/apt/* \
#            /var/lib/apt/lists/* \
#            /usr/share/man/?? \
#            /usr/share/man/??_*
