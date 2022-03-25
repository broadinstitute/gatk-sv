ARG VIRTUAL_ENV_IMAGE=sv-pipeline-children-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG GERT_DEPENDENCY_R_PKGS="openssl credentials askpass"
ARG SV_PIPELINE_R_PKGS="beeswarm devtools HardyWeinberg nloptr RColorBrewer vioplot zoo"

RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
        make automake pkg-config \
        libgit2-dev libssh2-1-dev libxml2-dev libssl-dev

# wget newer version of cmake because nloptr package crashes with the version available in ubuntu:18.04 (it needs -j -2)
RUN wget --quiet https://github.com/Kitware/CMake/releases/download/v3.21.5/cmake-3.21.5-linux-x86_64.tar.gz \
         -O /tmp/cmake.tar.gz && \
    tar -xf /tmp/cmake.tar.gz -C /tmp && \
    ln -s /tmp/cmake-3.21.5-linux-x86_64/bin/cmake /usr/local/bin

RUN install_R_packages.R $GERT_DEPENDENCY_R_PKGS
RUN install_deprecated_R_package.sh "https://cran.r-project.org/src/contrib/Archive/gert/gert_0.3.tar.gz"
# for some reason libxml2 is not discovered by R unless the path is specified manually
RUN export PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig/libxml-2.0.pc && \
    install_R_packages.R $SV_PIPELINE_R_PKGS

# Copy packages from build_virtual_env. Poofs the changes to the virtual environment without any extra cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
