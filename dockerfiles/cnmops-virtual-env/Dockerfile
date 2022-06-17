ARG VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG UBUNTU_RELEASE="22.04"

FROM $VIRTUAL_ENV_IMAGE

ARG RUN_DEPS="file libxml2"
ARG BUILD_DEPS="make cmake g++ gcc gfortran \
                libbz2-dev libblas-dev libicu-dev liblapack-dev liblzma-dev libpcre2-dev zlib1g-dev \
                libcurl4-openssl-dev libxml2-dev"
# update the list of required apt packages
RUN export NEW_PACKAGES=$(diff_of_lists.sh "$RUN_DEPS" $APT_REQUIRED_PACKAGES) && \
    printf " $NEW_PACKAGES" | fix_spaces.sh >> $APT_REQUIRED_PACKAGES

# Install R and bioconductor dependencies and packages
ARG DEBIAN_FRONTEND=noninteractive
RUN export APT_TRANSIENT_PACKAGES=$(diff_of_lists.sh "$BUILD_DEPS" $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends $BUILD_DEPS $RUN_DEPS && \
    install_bioconductor_packages.R cn.mops rtracklayer XML && \
    apt-get -qqy remove --purge $APT_TRANSIENT_PACKAGES && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
