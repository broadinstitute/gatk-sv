# This docker image contains the following
# list of tools and their dependencies:
# - ExpansionHunter (custom-built);
# - htslib, bcftools, and samtools;
# - python3.8, pip

FROM ubuntu:18.04

ARG HTSLIB_VER=1.15.1
ARG BCFTOOLS_VER=1.15.1
ARG SAMTOOLS_VER=1.15.1

RUN apt-get update && apt-get install --no-install-recommends -qqy \
    libcurl3 \
    tabix \
    git \
    make \
    gcc \
    g++ \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libncurses5-dev \
    libncursesw5-dev \
    python3.6-dev \
    python3-pip \
    python3-wheel \
    python3-venv \
    python3-setuptools \
    jq \
    wget \
    ca-certificates && \
    update-ca-certificates && \
    ln -s /usr/bin/python3.6 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip

RUN pip install --upgrade setuptools
# Combining the following two `RUN`s under a single `RUN` cause image build failure.
# RUN pip3 install -q -U wheel Cython==0.29.30 numpy==1.19.5 pandas==1.1.5
RUN pip install -q -U Cython==0.29.30 numpy==1.19.5
RUN pip install -q -U pandas==1.1.5

RUN mkdir eh && \
    wget -O eh/ExpansionHunter https://github.com/bw2/ExpansionHunter/raw/master/bin/linux/ExpansionHunter && \
    chmod 0755 /eh/ExpansionHunter
ENV PATH="/eh/:$PATH"

# Install htslib, bcftools, & samtools.
RUN wget -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/htslib-$HTSLIB_VER.tar.bz2 && \
    mkdir htslib && tar -xvjf htslib.tar.bz2 -C htslib --strip-components=1 && rm htslib.tar.bz2 && \
    cd htslib && ./configure && make && make install && cd .. && \
    wget -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/$BCFTOOLS_VER/bcftools-$BCFTOOLS_VER.tar.bz2 && \
    mkdir bcftools && tar -xvjf bcftools.tar.bz2 -C bcftools --strip-components=1 && rm bcftools.tar.bz2 && \
    cd bcftools && ./configure && make && make install && cd .. && \
    wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VER/samtools-$SAMTOOLS_VER.tar.bz2 && \
    mkdir samtools && tar -xvjf samtools.tar.bz2 -C samtools --strip-components=1 && rm samtools.tar.bz2 && \
    cd samtools && ./configure && make && make install && cd ..
ENV PATH="/bcftools:$PATH"
ENV PATH="/htslib:$PATH"
ENV PATH="/samtools:$PATH"

# Copy the scripts needed in the ExpansionHunter WDL.
COPY /src/str /opt/str

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/*
