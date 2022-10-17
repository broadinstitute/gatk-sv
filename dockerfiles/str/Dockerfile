# This docker image contains the following
# list of tools and their dependencies:
# - ExpansionHunter v5;
# - htslib, bcftools, and samtools;
# - python3.6, pip
# - numpy, and pandas

FROM ubuntu:22.10

ARG HTSLIB_VER=1.15.1
ARG BCFTOOLS_VER=1.15.1
ARG SAMTOOLS_VER=1.15.1
ARG EH_VER=5.0.0

RUN apt-get update && apt-get install --no-install-recommends -qqy \
    curl \
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
    build-essential \
    libssl-dev \
    libffi-dev \
    python3.10-dev \
    python3-pip \
    python3-wheel \
    python3-venv \
    python3-setuptools \
    jq \
    wget \
    ca-certificates && \
    update-ca-certificates && \
    ln -s /usr/bin/python3 /usr/bin/python

RUN pip install --upgrade setuptools
# Combining the following two `RUN`s under a single `RUN` cause image build failure.
# RUN pip3 install -q -U wheel Cython==0.29.30 numpy==1.19.5 pandas==1.1.5
RUN pip install -q numpy==1.23.4
RUN pip install -q Cython==0.29.32
RUN pip install -q pandas==1.5.0

RUN cd /opt/ && mkdir eh && cd eh/ && \
    wget -O eh.tar.gz https://github.com/Illumina/ExpansionHunter/releases/download/v$EH_VER/ExpansionHunter-v$EH_VER-linux_x86_64.tar.gz && \
    mkdir eh && tar -xf eh.tar.gz -C eh --strip-components=1 && \
    mv eh/bin/ExpansionHunter . && \
    rm eh.tar.gz && rm -r eh && cd ..
ENV PATH="/opt/eh/:$PATH"

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
ENV PATH="/opt/bcftools:$PATH"
ENV PATH="/opt/htslib:$PATH"
ENV PATH="/opt/samtools:$PATH"

# Copy the scripts needed in the ExpansionHunter WDL.
COPY /src/str /opt/str

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/*
