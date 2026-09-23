# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

############### stage 0: build samtools and bcftools from source
FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.9
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

ENV HOME=/root

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


# install GCloud SDK
ARG CLOUD_SDK_VERSION=239.0.0
ARG GCLOUD_SDK_EXTRA_COMPONENTS
ARG GCLOUD_SDK_TRANSIENT_DEP="curl gpg-agent gnupg python-pip python-setuptools "
ARG GCLOUD_SDK_DEPENDENCIES="lsb-release ${GCLOUD_SDK_TRANSIENT_DEP}"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${GCLOUD_SDK_DEPENDENCIES} && \
    pip install -q -U --no-cache-dir \
                crcmod \
                wheel && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends \
                 google-cloud-sdk=${CLOUD_SDK_VERSION}-0 \
                 $GCLOUD_SDK_EXTRA_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy purge ${GCLOUD_SDK_TRANSIENT_DEP} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*




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

# install vapor (optimized version, talkowski-lab/vapor) into its own Python 3.9 environment.
# Ubuntu 18.04 provides neither Python >= 3.9 nor the glibc that current pysam wheels need, so the
# interpreter and pysam come from conda-forge/bioconda through micromamba; the other libraries are
# installed at the versions pinned in vapor's requirements.txt (the versions it was tested with).
ARG VAPOR_REPO=https://github.com/talkowski-lab/vapor
ARG VAPOR_COMMIT=19e3812000f1553557eda45a37fa396b118b372c
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 ca-certificates curl && \
    rm -rf /var/lib/apt/lists/* && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /usr/local bin/micromamba && \
    micromamba create -y -p /opt/vapor-env -c conda-forge -c bioconda python=3.9 pip pysam=0.24.1 && \
    micromamba clean -a -y
RUN mkdir -p /opt/vapor && \
    curl -Ls ${VAPOR_REPO}/archive/${VAPOR_COMMIT}.tar.gz | tar -xz -C /opt/vapor --strip-components=1 && \
    /opt/vapor-env/bin/pip install --no-cache-dir -r /opt/vapor/requirements.txt && \
    /opt/vapor-env/bin/pip install --no-cache-dir --no-deps /opt/vapor && \
    ln -sf /opt/vapor-env/bin/vapor /usr/local/bin/vapor && \
    vapor bed --help | head -3

# parallelism comes from `vapor bed --threads`; one BLAS/OpenMP thread per worker process
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
