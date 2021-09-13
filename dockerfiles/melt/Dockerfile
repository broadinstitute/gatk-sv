FROM ubuntu:18.04

# Dependencies, including BOWTIE_2 and OpenJRE
ARG BOWTIE2_RELEASE=2.3.4.3
ARG BOWTIE2_DIR=/opt/bowtie2
ARG DEBIAN_FRONTEND=noninteractive
ENV JAVA_VERSION=8
ENV PATH="${BOWTIE2_DIR}:${PATH}"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 bzip2 \
                 ca-certificates \
                 libbz2-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 libtbb-dev \
                 libz-dev \
                 openjdk-${JAVA_VERSION}-jre-headless \
                 tabix \
                 unzip \
                 vcftools \
                 wget \
                 zlib1g-dev && \
    wget -q -P /tmp/ https://sourceforge.net/projects/bowtie-bio/files/bowtie2/${BOWTIE2_RELEASE}/bowtie2-${BOWTIE2_RELEASE}-linux-x86_64.zip && \
    cd /tmp/ && unzip -qq bowtie2-${BOWTIE2_RELEASE}-linux-x86_64.zip && \
    mkdir -p ${BOWTIE2_DIR} && mv bowtie2-${BOWTIE2_RELEASE}-linux-x86_64/* ${BOWTIE2_DIR} && \
    apt-get -qqy remove \
                 wget \
                 unzip && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Add melt jar
ARG MELT_RELEASE=2.0.5
ARG MELT_DIR="/MELT/"
RUN mkdir -p ${MELT_DIR}
WORKDIR ${MELT_DIR}
ADD MELTv${MELT_RELEASE}_patch.tar.gz .
RUN echo "export PATH=${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH" > /etc/profile.d/MELT.sh
ENV PATH="${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH"

# Add master script
COPY dockerfiles/melt/run_MELT_${MELT_RELEASE}.sh .
RUN chmod +x run_MELT_${MELT_RELEASE}.sh
