ARG UBUNTU_RELEASE="22.04"
# available R releases listed here:  https://cdn.rstudio.com/r/versions.json
ARG R_RELEASE_VERSION="4.5.0"
ARG R_INSTALL_PATH=/opt/R
ARG APT_REQUIRED_PACKAGES="/opt/apt-required-packages.list"

######################## install R
FROM ubuntu:$UBUNTU_RELEASE as build_r
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR
ARG CONDA_RELEASE
ARG VIRTUAL_ENV_PATH

# copy apt helper script and update path
RUN mkdir -p /opt/bin
COPY dockerfiles/sv-base-virtual-env/fix_spaces.sh /opt/bin
ENV PATH=/opt/bin:$PATH

# get the list of all packages needed to support R (and the virtual env in general)
ARG BUILD_DEPS="make cmake g++ gcc gfortran \
                libbz2-dev libopenblas-dev libicu-dev liblapack-dev liblzma-dev libpcre2-dev zlib1g-dev"
ARG RUN_DEPS="wget curl apt-transport-https ca-certificates unzip zip ucf \
              libbz2-1.0 libopenblas0 libicu70 liblapack3 liblzma5 libpcre2-8-0 libpcre2-16-0 libpcre2-32-0 \
              libpcre2-posix3 zlib1g libc6 libcairo2 libcurl4 libglib2.0-0 libgomp1 libjpeg8 libpango-1.0-0 \
              libpangocairo-1.0-0 libpaper-utils libpng16-16 libreadline8 libtcl8.6 libtiff5 libtk8.6 libx11-6 libxt6 \
              libdeflate-dev libssl-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libwebp-dev"
ARG APT_REQUIRED_PACKAGES
RUN printf "$RUN_DEPS" | fix_spaces.sh > $APT_REQUIRED_PACKAGES

# use apt to install transient dependencies for installing R
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends $BUILD_DEPS $RUN_DEPS

# download the .deb file for installing R (heavily modified from instructions at https://github.com/rstudio/r-builds)
ARG UBUNTU_RELEASE
ARG R_RELEASE_VERSION
ARG R_DEB="r-${R_RELEASE_VERSION}_1_amd64.deb"
RUN curl --connect-timeout 10 -O https://cdn.rstudio.com/r/ubuntu-$(echo $UBUNTU_RELEASE | sed 's/\.//')/pkgs/$R_DEB

# install R, and symlink to a standard location
ARG R_INSTALL_PATH
ARG R_INSTALL_BIN=$R_INSTALL_PATH/bin
RUN dpkg -i $R_DEB && \
    ln -s /opt/R/$R_RELEASE_VERSION/* $R_INSTALL_PATH/

# copy R helper scripts and update path
RUN mkdir -p $R_INSTALL_BIN
COPY dockerfiles/sv-base-virtual-env/install_R_packages.R $R_INSTALL_BIN/
COPY dockerfiles/sv-base-virtual-env/install_bioconductor_packages.R $R_INSTALL_BIN/
COPY dockerfiles/sv-base-virtual-env/install_deprecated_R_package.sh $R_INSTALL_BIN/
COPY dockerfiles/sv-base-virtual-env/rm_unneeded_r_library_files.sh $R_INSTALL_BIN/
COPY dockerfiles/sv-base-virtual-env/diff_of_lists.sh $R_INSTALL_BIN/
COPY dockerfiles/sv-base-virtual-env/union_lists.sh $R_INSTALL_BIN/
ENV PATH=$R_INSTALL_BIN:/opt/bin:$PATH

# update the installed R packages
RUN Rscript -e "update.packages(contriburl=contrib.url(repos=c('http://lib.stat.cmu.edu/R/CRAN/', 'https://cran.rstudio.com')), clean=TRUE, quiet=TRUE, ask=FALSE)"
# install the required R packages for this stage
ARG SV_BASE_R_PKGS="optparse BiocManager"
RUN install_R_packages.R $SV_BASE_R_PKGS


########################## Copy R install to final image
# Effectively removes everything but the minimal environment and flattens docker layers in the final image.
FROM ubuntu:$UBUNTU_RELEASE

# copy R install
COPY --from=build_r /opt /opt

# update paths and ENV variables
ARG R_INSTALL_PATH
ARG R_INSTALL_BIN=$R_INSTALL_PATH/bin
ENV PATH=$R_INSTALL_BIN:/opt/bin:$PATH
ARG R_RELEASE_VERSION
ENV R_RELEASE_VERSION=$R_RELEASE_VERSION
ARG APT_REQUIRED_PACKAGES
ENV APT_REQUIRED_PACKAGES=$APT_REQUIRED_PACKAGES

# install needed support packages
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends $(fix_spaces.sh $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
