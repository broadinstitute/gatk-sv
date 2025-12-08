ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
ARG WHAM_IMAGE=wham:latest

FROM $WHAM_IMAGE as wham_layer
FROM $SV_PIPELINE_IMAGE

ENV DEBIAN_FRONTEND noninteractive

# apt-get update and install global requirements
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autoconf \
        autogen \
        build-essential \
        curl \
        git \
        libbz2-dev \
        libcurl4-openssl-dev \
        liblzma-dev \
        libncurses5-dev \
        libnss-sss \
        libssl-dev \
        libxml2-dev \
        ncbi-blast+ \
        r-base \
        r-bioc-biostrings \
        r-bioc-rsamtools \
        r-cran-biocmanager \
        r-cran-devtools \
        r-cran-stringr \
        r-cran-optparse \
        wget \
        zlib1g-dev \
        jq \
        # needed for manta
        python2.7

# apt-get clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install global r requirements
RUN Rscript -e "install.packages('htmlwidgets', repos='http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('shiny', repos='http://cran.us.r-project.org')"
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "library(devtools); install_github('mhahsler/rBLAST')"

# install bcftools and htslib
ARG HTSLIB_VERSION="1.18"
RUN mkdir -p /opt && cd /opt && \
    wget -q https://github.com/samtools/bcftools/releases/download/${HTSLIB_VERSION}/bcftools-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf bcftools-${HTSLIB_VERSION}.tar.bz2 && \
    cd bcftools-${HTSLIB_VERSION} && \
    ./configure --quiet && \
    make -s all && \
    make -s install install-htslib && \
    cd / && rm -r /opt/bcftools-${HTSLIB_VERSION} /opt/bcftools-${HTSLIB_VERSION}.tar.bz2

ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# install manta
ARG MANTA_RELEASE=1.6.0
ARG MANTA_INSTALL_DIR=/usr/local/bin/manta/
RUN wget -q https://github.com/Illumina/manta/releases/download/v${MANTA_RELEASE}/manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    tar -xjf manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    mv ./manta-${MANTA_RELEASE}.centos6_x86_64 ${MANTA_INSTALL_DIR} && \
    ln -s $(which python2.7) /usr/bin/python2 && \
    python2 ${MANTA_INSTALL_DIR}/bin/runMantaWorkflowDemo.py && \
    rm -rf MantaDemoAnalysis && \
    rm -rf ${MANTA_INSTALL_DIR}share/demo

# install scramble
ARG SCRAMBLE_COMMIT="56b5ae849d16ec1fc83ea1426b0ffc356ee6d99c"
RUN mkdir /app && cd /app \
    && git clone https://github.com/mwalker174/scramble-gatk-sv.git \
    && cd scramble-gatk-sv \
    && git checkout ${SCRAMBLE_COMMIT} \
    && cd cluster_identifier/src \
    && make \
    && ln -s /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier /usr/local/bin

# test
RUN Rscript --vanilla /app/scramble-gatk-sv/cluster_analysis/bin/SCRAMble.R --help
RUN /app/scramble-gatk-sv/cluster_identifier/src/build/cluster_identifier -v

# define default command
CMD ["Rscript"]


# Whamg Docker
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 gawk \
                 libgomp1 \
                 tabix \
                 zlib1g-dev && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean
COPY --from=wham_layer /bin/whamg /bin/


# cnmops
# Note that the following installation is different from the installation in cnmops docker
# since that method does not compile on this docker.
RUN mkdir -p /opt/R/4.1.3/lib/R/doc/html && \
    R -e 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/focal/latest")); BiocManager::install(c("cn.mops", "rtracklayer", "XML"), update=FALSE, ask=FALSE)'


# Install GATK
# note that gatk is installed in sv-base docker, which this docker image inherits.
# however, the base images are currently being updated: https://github.com/broadinstitute/gatk-sv/pull/816
# so as a temporary work-around we're re-installing gatk in this image using code in a specific branch.
# setting SHELL to bash in the following is needed for the `source` command in the following RUN,
# since the default sh does not support `source`.
SHELL ["/bin/bash", "-c"]
RUN apt-get -qqy install --no-upgrade --no-install-recommends git git-lfs openjdk-17-jdk openjdk-17-jre-headless libgomp1 && \
    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && \
    git checkout mw_gatk_sv && \
    git checkout 1eac9db89b65f9493d2a91539e2b9274159a5b5c && \
    mkdir -p /opt/gatk_miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh -O /opt/gatk_miniconda3/miniconda.sh && \
    bash /opt/gatk_miniconda3/miniconda.sh -p /opt/gatk_miniconda3 -b -u && \
    rm /opt/gatk_miniconda3/miniconda.sh && \
    /opt/gatk_miniconda3/bin/conda config --append envs_dirs /opt/conda/envs && \
    source /opt/gatk_miniconda3/bin/activate && \
    conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    ./gradlew localDevCondaEnv && \
    ./gradlew localJar && \
    cp ./build/libs/gatk.jar /opt/gatk.jar && \
    cd .. && \
    rm -rf gatk/ && \
    ln -s /opt/gatk_miniconda3/bin/conda /usr/local/bin/conda


WORKDIR /