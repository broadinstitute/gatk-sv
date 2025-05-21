ARG SV_PIPELINE_IMAGE=sv-pipeline:latest

FROM us.gcr.io/broad-dsde-methods/gatk-sv/wham:2024-10-25-v0.29-beta-5ea22a52 as wham_layer

# FROM $SV_PIPELINE_IMAGE
FROM us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-02-11-v1.0.2-hotfix-22bf77e0

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

# FROM us.gcr.io/broad-dsde-methods/gatk-sv/manta:2023-09-14-v0.28.3-beta-3f22f94d as manta_layer

# FROM us.gcr.io/broad-dsde-methods/gatk-sv/gatk:2024-12-05-4.6.1.0-6-gfc248dfc1-NIGHTLY-SNAPSHOT as gatk_layer

WORKDIR /