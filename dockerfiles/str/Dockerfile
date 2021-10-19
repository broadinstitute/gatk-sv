# This docker image contains the following
# list of tools and their dependencies:
# - GangSTR
# - TRTools
# - ExpansionHunter

FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG CONDA_RELEASE="4.6.14"
ARG CONDA_INSTALL_DIR="/opt/conda"
ARG CONDA_BIN=${CONDA_INSTALL_DIR}/bin
ARG CONDA_CMD=${CONDA_BIN}/conda
ENV PATH=${CONDA_BIN}:$PATH
ENV EH_VERSION=v4.0.2

RUN apt-get update && apt-get install --no-install-recommends -qqy \
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
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev

# Install and configure conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_RELEASE}-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_INSTALL_DIR} && \
    ${CONDA_CMD} init bash && \
    bash -c "source ${HOME}/.bashrc"

RUN pip install Cython==0.29.24 numpy==1.21.2 pybedtools==0.8.2 scipy==1.7.1

# Install samtools (needed to index reference fasta files)
RUN wget -O samtools-1.13.tar.bz2 https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xjf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure --without-curses && make && make install && \
    cd ..

# Install bedtools (needed for DumpSTR)
## Option 1: install from source
RUN wget -O bedtools-2.27.1.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
RUN tar -xzvf bedtools-2.27.1.tar.gz
WORKDIR bedtools2
RUN make && make install
WORKDIR ..
## Option 2: install from apt
#RUN apt-get install bedtools

# Download and install TRTools
RUN git clone https://github.com/gymreklab/TRTools && \
    cd TRTools && \
    git checkout tags/v4.0.0 -b v4.0.0-branch && \
    pip install -r requirements.txt && \
    python3 setup.py install && \
    cd ..

# Install ExpansionHunter
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/${EH_VERSION}/ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    tar xzf ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    rm ExpansionHunter-${EH_VERSION}-linux_x86_64.tar.gz && \
    mv /ExpansionHunter-${EH_VERSION}-linux_x86_64 /ExpansionHunter
ENV PATH="/ExpansionHunter/bin/:$PATH"

# Install GangSTR
RUN conda install -c bioconda -c conda-forge gangstr

# Clean up
RUN rm -rf /tmp/* \
    /var/tmp/* \
    /var/cache/apt/* \
    /var/lib/apt/lists/* \
