############### stage 0: build samtools from source
FROM ubuntu:20.04

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

############### stage 1: install Manta
FROM ubuntu:20.04

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/

# install necessary packages while building, and tools necessary for running, then clean up; all in one layer
# note that we explicitly run the demo/test-run distributed with the binary, then delete all demo results and data
ARG MANTA_RELEASE=1.6.0
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
