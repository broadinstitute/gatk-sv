# source Image
FROM ubuntu:22.04

# set noninterative mode
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
        zlib1g-dev

# apt-get clean and remove cached source lists
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install global r requirements
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "options(repos = 'https://CRAN.R-project.org'); library(devtools); install_github('mhahsler/rBLAST')"

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
