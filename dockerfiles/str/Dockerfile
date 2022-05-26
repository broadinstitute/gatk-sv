# This docker image contains the following
# list of tools and their dependencies:
# - ExpansionHunter (custom-built)
# - htslib, bcftools, and samtools.

FROM ubuntu:18.04

ARG HTSLIB_VER=1.15.1
ARG BCFTOOLS_VER=1.15.1

RUN apt-get update && apt-get install --no-install-recommends -qqy \
    libcurl3 \
    tabix \
    git \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libncurses5-dev \
    libncursesw5-dev \
    jq \
    wget \
    ca-certificates \
    && update-ca-certificates

RUN mkdir eh && \
    wget -O eh/ExpansionHunter https://github.com/bw2/ExpansionHunter/raw/master/bin/linux/ExpansionHunter && \
    chmod 0755 /eh/ExpansionHunter
ENV PATH="/eh/:$PATH"

RUN wget -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/htslib-$HTSLIB_VER.tar.bz2 && \
    mkdir htslib && tar -xvjf htslib.tar.bz2 -C htslib --strip-components=1 && rm htslib.tar.bz2 && \
    cd htslib && ./configure && make && make install && cd .. && \
    wget -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2 && \
    mkdir bcftools && tar -xvjf bcftools.tar.bz2 -C bcftools --strip-components=1 && rm bcftools.tar.bz2 && \
    cd bcftools && make && make install && cd .. && \
    wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar -xvjf samtools-1.15.1.tar.bz2 && \
    rm samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure && \
    make && \
    make install

ENV PATH="/bcftools-1.15.1:$PATH"
ENV PATH="/htslib-1.15.1:$PATH"
ENV PATH="/samtools-1.15.1:$PATH"

COPY /src/str /opt/str

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/*
