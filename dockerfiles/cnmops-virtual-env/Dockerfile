ARG VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG UBUNTU_RELEASE="18.04"

FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

# Install R and bioconductor dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake libcurl4-openssl-dev \
                 file g++ gfortran liblapack-dev libopenblas-dev libxml2-dev && \
    install_deprecated_R_package.sh https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz && \
    install_bioconductor_packages.R cn.mops rtracklayer

# Copy packages from build_virtual_env. Poofs the changes to the virtual environment without any extra cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
