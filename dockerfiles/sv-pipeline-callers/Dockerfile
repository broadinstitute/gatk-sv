# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_IMAGE=gatksv/sv-pipeline-base:latest
FROM ${SV_PIPELINE_BASE_IMAGE}

#########################################################
### MANTA
#########################################################

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.9
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
    cd samtools-${SAMTOOLS_VERSION} && ./configure CFLAGS="-static" --without-curses && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help

# install necessary packages while building, and tools necessary for running, then clean up; all in one layer
# note that we explicitly run the demo/test-run distributed with the binary, then delete all demo results and data
ARG MANTA_RELEASE=1.5.0
ARG MANTA_INSTALL_DIR=/usr/local/bin/manta/
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 ca-certificates \
                 python2.7 \
                 tabix \
                 wget \
                 zlib1g-dev && \
    cd /opt && \
    wget -q https://github.com/Illumina/manta/releases/download/v${MANTA_RELEASE}/manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    tar -xjf manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    mv /opt/manta-${MANTA_RELEASE}.centos6_x86_64 ${MANTA_INSTALL_DIR} && \
    ln -s $(which python2.7) /usr/bin/python2 && \
    python2 ${MANTA_INSTALL_DIR}/bin/runMantaWorkflowDemo.py && \
    rm -rf MantaDemoAnalysis && \
    rm -rf ${MANTA_INSTALL_DIR}share/demo && \
    apt-get -qqy remove \
                 wget \
                 ca-certificates && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 \
           /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help


#########################################################
### MELT
#########################################################

# Dependencies, including BOWTIE_2 and OpenJRE
ARG BOWTIE2_RELEASE=2.3.4.3
ARG BOWTIE2_DIR=/opt/bowtie2
ARG DEBIAN_FRONTEND=noninteractive
ENV JAVA_VERSION=8
ENV PATH="${BOWTIE2_DIR}:${PATH}"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 bzip2 \
                 ca-certificates \
                 libbz2-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 libtbb-dev \
                 libz-dev \
                 openjdk-${JAVA_VERSION}-jre-headless \
                 tabix \
                 unzip \
                 vcftools \
                 wget \
                 zlib1g-dev && \
   wget --no-check-certificate -q -P /tmp/ https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_RELEASE}/bowtie2-${BOWTIE2_RELEASE}-linux-x86_64.zip && \
   cd /tmp/ && unzip -qq bowtie2-${BOWTIE2_RELEASE}-linux-x86_64.zip && \
   mkdir -p ${BOWTIE2_DIR} && mv bowtie2-${BOWTIE2_RELEASE}-linux-x86_64/* ${BOWTIE2_DIR} && \
   apt-get -qqy remove \
                 wget \
                 unzip && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Add melt jar
ARG MELT_RELEASE=2.0.5
ARG MELT_DIR="/MELT/"
RUN mkdir -p ${MELT_DIR}
WORKDIR ${MELT_DIR}
ADD MELTv${MELT_RELEASE}_patch.tar.gz .
RUN echo "export PATH=${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH" > /etc/profile.d/MELT.sh
ENV PATH="${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH"

# Add master script
COPY run_MELT_${MELT_RELEASE}.sh .
RUN chmod +x run_MELT_${MELT_RELEASE}.sh

#########################################################
### WHAM
#########################################################

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
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

ADD whamg /bin/


#########################################################
### PICARD
#########################################################

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends wget && \
    mkdir -p /opt/picard && \
    cd /opt/picard && \
    wget -q --no-check-certificate https://github.com/broadinstitute/picard/releases/download/2.25.0/picard.jar && \
    apt-get -qqy purge wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*


#########################################################
### R
#########################################################

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends r-base-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
