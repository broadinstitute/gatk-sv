# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

############### stage 0: build samtools and bcftools from source
FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.11
ARG BCFTOOLS_VERSION=1.9
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help

############### stage 1: other commonly used bioinformatics utilities 
FROM ubuntu:18.04

# Set environment variables.
ENV HOME=/root
ENV R_HOME=/usr


# Define default command.
CMD ["bash"]

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/

#### Basic utilities
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bedtools \
                 curl \
                 gawk \
                 less \
                 tabix \
                 vcftools \
                 wget \
                 zlib1g-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help


# GATK needed for PrintReads
ARG GATK_TRANSIENT_DEP="git git-lfs"
ARG GATK_DEP="openjdk-8-jdk ${GATK_TRANSIENT_DEP}"
ARG GATK_COMMIT="a33bf19dd3188af0af1bd17bce015eb20ba73227"
ENV GATK_JAR="/opt/gatk.jar"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${GATK_DEP} && \
    cd / && git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && git lfs install && \
    git checkout ${GATK_COMMIT} && ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR} && \
    cd / && \
    apt-get -qqy purge ${GATK_TRANSIENT_DEP} && \
    apt-get -qqy clean && \
    rm -rf /gatk \
           /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    java -jar ${GATK_JAR}

# Define default command.
CMD ["bash"]

#### Basic image utilities

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get update && apt-get install -y r-base
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y \
    bcftools \
    bedtools \
    openjdk-8-jdk \
    xvfb \
    zip \
    python3-pip \
    python3-dev \
    curl \
    git \
    samtools \
    tabix \
    vcftools \
    wget \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get -y clean && \
    apt-get -y autoremove && \
    apt-get  -y autoclean

# install vapor
RUN pip3 install --no-cache-dir Cython && \
    cd /opt && \
    git clone https://github.com/mills-lab/vapor.git && \
    cd /opt/vapor/ && \
    git checkout 5e322b5 && \
    pip3 install -e .



