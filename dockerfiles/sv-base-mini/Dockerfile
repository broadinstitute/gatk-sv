# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.
ARG UBUNTU_RELEASE="22.04"
ARG HTSLIB_VERSION="1.15.1"
ARG BEDTOOLS_VERSION="2.31.0"
ARG VCFTOOLS_VERSION="0.1.16"
ARG BWA_COMMIT="139f68fc4c3747813783a488aef2adc86626b01b"

ARG APT_REQUIRED_PACKAGES="/opt/apt-required-packages.list"

############### stage 0: download and build tools
FROM ubuntu:$UBUNTU_RELEASE as build

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
        ca-certificates autoconf automake bzip2 g++ git make wget pkgconf python2 \
        libssl-dev libbz2-dev libcurl4-openssl-dev liblzma-dev libncurses-dev zlib1g-dev libdeflate-dev

# install samtools
ARG HTSLIB_VERSION
RUN wget -q https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf samtools-${HTSLIB_VERSION}.tar.bz2 && \
    cd samtools-${HTSLIB_VERSION} && \
    ./configure --quiet --enable-libcurl --prefix=/opt/samtools && \
    make -s -j $(nproc --all) all all-htslib && \
    make -s install install-htslib
# install bcftools
RUN wget -q https://github.com/samtools/bcftools/releases/download/${HTSLIB_VERSION}/bcftools-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf bcftools-${HTSLIB_VERSION}.tar.bz2 && \
    cd bcftools-${HTSLIB_VERSION} && \
    ./configure --quiet --prefix=/opt/bcftools && \
    make -s -j $(nproc --all) && \
    make -s install
# install vcftools
ARG VCFTOOLS_VERSION
RUN wget -q https://github.com/vcftools/vcftools/releases/download/v$VCFTOOLS_VERSION/vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    tar zxf vcftools-${VCFTOOLS_VERSION}.tar.gz && \
    cd vcftools-${VCFTOOLS_VERSION} && \
    export PERL5LIB=$(pwd)/src/perl/ && \
    ./configure --quiet --prefix=/opt/vcftools && \
    make -s -j $(nproc --all) && \
    make -s install

# install bedtools
#   during install, bedtools needs python2 to be callable as "python"
RUN ln -s /usr/bin/python2 /usr/bin/python
ARG BEDTOOLS_VERSION
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VERSION/bedtools.static && \
    mkdir -p /opt/bedtools/bin && \
    mv bedtools.static /opt/bedtools/bin/bedtools && \
    chmod a+x /opt/bedtools/bin/bedtools

# install bwa
# must do from source because of compiler error in latest release (see https://github.com/lh3/bwa/issues/387)
ARG BWA_COMMIT
RUN cd /opt && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    git checkout $BWA_COMMIT && \
    make -s && \
    cd .. && \
    mkdir -p /opt/bin && \
    mv /opt/bwa/bwa /opt/bin/ && \
    rm -r bwa
ENV PATH=/opt/bin:$PATH

############### stage 1: copy tools and install needed non-dev libraries
FROM ubuntu:$UBUNTU_RELEASE

ENV HOME=/root

# Define default command.
CMD ["bash"]

# copy utility file for cleaning up $APT_REQUIRED_PACKAGES
RUN mkdir -p /opt/bin
COPY dockerfiles/sv-base-mini/fix_spaces.sh /opt/bin/
# copy the results of the build stage, and update path
COPY --from=build /opt /opt
ENV PATH=/opt/samtools/bin:/opt/bcftools/bin:/opt/vcftools/bin:/opt/bedtools/bin:/opt/bin:$PATH
# note: need to specify path to perl libraries or tools from vcftools won't work
RUN LIB_LOC=$(find /opt/vcftools/share/perl -mindepth 1 -type d) && \
    ln -s $LIB_LOC /opt/vcftools/share/perl_lib
ENV PERL5LIB=/opt/vcftools/share/perl_lib

# establish file storing required apt packages (so subsequent stages know which packages can be removed and which can't)
ARG APT_REQUIRED_PACKAGES
ENV APT_REQUIRED_PACKAGES=$APT_REQUIRED_PACKAGES
ARG REQUIRED_DEPS="apt-transport-https ca-certificates curl wget gawk less \
                   libssl3 libbz2-1.0 liblzma5 libcurl4 libncurses6 zlib1g libdeflate0 libdata-dumper-compact-perl"
RUN printf "$REQUIRED_DEPS" | fix_spaces.sh > $APT_REQUIRED_PACKAGES

# upgrade OS & original packages, and install required packages
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends $REQUIRED_DEPS && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Make the prompt color. Mildly convenient for debugging purposes
RUN sed -i 's/^#force_color_prompt=yes$/force_color_prompt=yes/' ~/.bashrc

# Show the installed software works
RUN samtools --version
RUN tabix --version
RUN bcftools --version
RUN bedtools --version
RUN vcftools --version
RUN which bwa
