ARG SVBASE_IMAGE=sv-base:latest
FROM $SVBASE_IMAGE

###############
# MELT
###############

# Dependencies, including BOWTIE_2 and OpenJRE
ARG BOWTIE2_RELEASE=2.3.4.3
ARG BOWTIE2_DIR=/opt/bowtie2
ARG DEBIAN_FRONTEND=noninteractive
ENV JAVA_VERSION=8
ENV PATH="${BOWTIE2_DIR}:${PATH}"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bzip2 \
                 ca-certificates \
                 libbz2-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 libtbb-dev \
                 libz-dev \
                 openjdk-${JAVA_VERSION}-jre-headless \
                 unzip \
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
ADD dockerfiles/melt/MELTv${MELT_RELEASE}_patch.tar.gz .
RUN echo "export PATH=${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH" > /etc/profile.d/MELT.sh
ENV PATH="${MELT_DIR}/MELTv${MELT_RELEASE}_patch/:$PATH"

# Add master script
COPY dockerfiles/melt/run_MELT_${MELT_RELEASE}.sh .
RUN chmod +x run_MELT_${MELT_RELEASE}.sh


###############
# Manta
###############

# install necessary packages while building, and tools necessary for running, then clean up; all in one layer
# note that we explicitly run the demo/test-run distributed with the binary, then delete all demo results and data
ARG MANTA_RELEASE=1.6.0
ARG MANTA_INSTALL_DIR=/usr/local/bin/manta/
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 bcftools \
                 ca-certificates \
                 python2.7 \
                 tabix \
                 wget \
                 zlib1g-dev && \
    cd /opt && \
    wget -q https://github.com/Illumina/manta/releases/download/v${MANTA_RELEASE}/manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    tar -xjf manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 && \
    mv /opt/manta-${MANTA_RELEASE}.centos6_x86_64 ${MANTA_INSTALL_DIR} && \
    ln -s $(which python2.7) /usr/bin/python2 && \
    python2 ${MANTA_INSTALL_DIR}/bin/runMantaWorkflowDemo.py && \
    rm -rf MantaDemoAnalysis && \
    rm -rf ${MANTA_INSTALL_DIR}share/demo && \
    apt-get -qqy remove \
                 wget \
                 ca-certificates && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf manta-${MANTA_RELEASE}.centos6_x86_64.tar.bz2 \
           /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help


###############
# Wham
###############

COPY dockerfiles/wham/whamg /bin/


###############
# GATK
###############

ARG GATK_COMMIT="c9bf9411c4531ebadb6c8c24d80e6c170e8b9691"
ARG GATK_JAR="/opt/gatk.jar"

ARG GATK_BUILD_DEP="git git-lfs openjdk-17-jdk"
ARG GATK_RUN_DEP="openjdk-17-jre-headless libgomp1"
ARG DEBIAN_FRONTEND=noninteractive
# install build deps
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 $GATK_BUILD_DEP $GATK_RUN_DEP
# clone requested GATK_COMMIT
RUN cd /opt && GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && \
    git checkout ${GATK_COMMIT}
# build GATK_JAR
RUN cd gatk && \
    ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR}
