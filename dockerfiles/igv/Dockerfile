# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

# Start with linux base
FROM ubuntu:18.04

# Set environment variables.
ENV HOME=/root



#### Basic image utilities
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

ARG CLOUD_SDK_VERSION=272.0.0
ENV CLOUD_SDK_VERSION=$CLOUD_SDK_VERSION


ARG INSTALL_COMPONENTS
RUN apt-get update -qqy && apt-get install -qqy \
        curl \
        gcc \
        python3-dev \
        python3-setuptools \
        apt-transport-https \
        lsb-release \
        openssh-client \
        git \
        gnupg \
    && pip3 install -U crcmod && \
    pip3 install numpy && \
    pip3 install image  && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version


RUN wget http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip && \
    unzip IGV_2.4.14.zip


VOLUME ["/root/.config"]
COPY MakeRDtest.py /src/MakeRDtest.py
COPY makeigvpesr_cram.py  /src/makeigvpesr_cram.py
COPY makeigvsplit_cram.py /src/makeigvsplit_cram.py
COPY makeigv_cram.py      /src/makeigv_cram.py
COPY makeigvpesr_trio.py  /src/makeigvpesr_trio.py
COPY prefs.properties     /root/igv/prefs.properties
COPY arial.ttf /usr/src/app/arial.ttf

# Define default command.
CMD ["bash"]