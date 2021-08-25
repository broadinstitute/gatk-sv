# This docker image contains the following
# list of tools and their dependencies:
# - GangSTR
# - TRTools
# - ExpansionHunter

FROM ubuntu:20.04

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install --no-install-recommends -qqy \
    python3-dev \
    python3-pip \
    python \
    python-dev \
    awscli \
    build-essential \
    git \
    libbz2-dev \
    liblzma-dev \
    make \
    pkg-config \
    wget \
    unzip \
    zlib1g-dev

RUN pip3 install pybedtools==0.8.2 pyvcf==0.6.8 scipy==1.7.1 numpy==1.21.1

# Install samtools (needed to index reference fasta files)
RUN wget -O samtools-1.9.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar -xjf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure --without-curses && make && make install \
    && cd ..

# Install bedtools (needed for DumpSTR)
## Option 1: install from source
RUN wget -O bedtools-2.27.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
RUN tar -xzvf bedtools-2.27.1.tar.gz
WORKDIR bedtools2
RUN make && make install
WORKDIR ..
## Option 2: install from apt
#RUN apt-get install bedtools

# Download, compile, and install GangSTR
RUN wget -O GangSTR-2.4.tar.gz https://github.com/gymreklab/GangSTR/releases/download/v2.4/GangSTR-2.4.tar.gz \
    && tar -xzvf GangSTR-2.4.tar.gz \
    && cd GangSTR-2.4 \
    && ./install-gangstr.sh \
    && ldconfig \
    && cd ..

# Download and install TRTools
RUN git clone https://github.com/gymreklab/TRTools \
    && cd TRTools \
    && python3 setup.py install \
    && cd ..

ENV EH_VERSION=v4.0.2
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
    && tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
    && rm ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz \
    && mv /ExpansionHunter-${EH_VERSION}-linux_x86_64 /ExpansionHunter
ENV PATH="/ExpansionHunter/bin/:$PATH"
