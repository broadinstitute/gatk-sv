# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

############### stage 0: build samtools and bcftools from source
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
# htslib 1.9 cannot fetch remote (gs://, https://) indexes with Ubuntu 22.04's libcurl; 1.24 can,
# and supports samtools view -X for an explicitly given index
ARG SAMTOOLS_VERSION=1.24
ARG BCFTOOLS_VERSION=1.24
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
FROM ubuntu:22.04

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
# google-cloud-sdk 239 is no longer published; the package is now google-cloud-cli (latest only)
ARG GCLOUD_SDK_EXTRA_COMPONENTS
ARG GCLOUD_SDK_TRANSIENT_DEP="curl gpg-agent gnupg python3-pip python3-setuptools "
ARG GCLOUD_SDK_DEPENDENCIES="lsb-release ${GCLOUD_SDK_TRANSIENT_DEP}"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ${GCLOUD_SDK_DEPENDENCIES} && \
    pip3 install -q -U --no-cache-dir \
                crcmod \
                wheel && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends \
                 google-cloud-cli \
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

# install vapor (optimized, talkowski-lab/vapor at VAPOR_COMMIT) into its own Python 3.6 environment.
# Scores depend on the numerical libraries, so they are pinned (requirements.txt) to the versions of the
# original Ubuntu 18.04 / Python 3.6 image; Ubuntu 22.04 has no Python 3.6, so micromamba provides it.
# The scoring module is compiled with Cython 0.29 like the original.
ARG VAPOR_REPO=https://github.com/talkowski-lab/vapor.git
ARG VAPOR_COMMIT=d7930417ebbb4ec1709a9dd5e4aa6fd9b0bc68a2
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 ca-certificates curl gcc git && \
    rm -rf /var/lib/apt/lists/* && \
    MAMBA_PLATFORM=$(case $(uname -m) in aarch64) echo linux-aarch64 ;; *) echo linux-64 ;; esac) && \
    curl -Ls https://micro.mamba.pm/api/micromamba/${MAMBA_PLATFORM}/latest | tar -xj -C /usr/local bin/micromamba && \
    micromamba create -y -p /opt/vapor-env -c conda-forge python=3.6 pip && \
    micromamba clean -a -y && \
    git clone ${VAPOR_REPO} /opt/vapor && \
    cd /opt/vapor && git checkout ${VAPOR_COMMIT} && \
    /opt/vapor-env/bin/pip install --no-cache-dir -r requirements.txt && \
    CC=gcc LDSHARED="gcc -pthread -shared" /opt/vapor-env/bin/pip install --no-cache-dir --no-deps --no-build-isolation . && \
    ln -sf /opt/vapor-env/bin/vapor /usr/local/bin/vapor && \
    cd / && /opt/vapor-env/bin/python -c "import vapor_vali.Simple_function as m; print(m.__file__)"

# parallelism comes from `vapor bed --threads`; one BLAS/OpenMP thread per worker process
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
