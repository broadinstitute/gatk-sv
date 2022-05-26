# This docker image contains the following
# list of tools and their dependencies:
# - ExpansionHunter
# - htslib, bcftools, and samtools.

FROM ubuntu:18.04

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

COPY --chmod=0755 /dockerfiles/str/ExpansionHunter /ExpansionHunter/
ENV PATH="/ExpansionHunter/:$PATH"


RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 && \
    tar -xvjf htslib-1.15.1.tar.bz2 && \
    rm htslib-1.15.1.tar.bz2 && \
    cd htslib-1.15.1 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2 && \
    tar -xvjf bcftools-1.15.1.tar.bz2 && \
    rm bcftools-1.15.1.tar.bz2 && \
    cd bcftools-1.15.1 && \
    make && \
    make install && \
    cd .. && \
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
