# GATK SV Pipeline RdTest virtual env dockerfile
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-children-virtual-env:latest
FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG DATA_TABLE_URL="https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.12.8.tar.gz"
ARG MNORMT_URL="https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
ARG DAAG_URL="https://cran.r-project.org/src/contrib/Archive/DAAG/DAAG_1.22.tar.gz"

ARG SV_PIPELINE_BIOC_PKGS="multtest"
ARG SV_PIPELINE_R_PKGS="BSDA caret fpc hash metap perm plyr pwr reshape ROCR zoo"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake pkg-config \
                 fftw3 fftw3-dev libssh2-1-dev libssl-dev libjpeg-dev libxml2-dev

RUN install_deprecated_R_package.sh $DATA_TABLE_URL
RUN install_R_packages.R "latticeExtra"
RUN install_bioconductor_packages.R $SV_PIPELINE_BIOC_PKGS
RUN install_deprecated_R_package.sh $MNORMT_URL
RUN install_deprecated_R_package.sh $DAAG_URL
# for some reason libxml2 is not discovered by R unless the path is specified manually
RUN export PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig/libxml-2.0.pc && \
    install_R_packages.R $SV_PIPELINE_R_PKGS

# Copy packages from build_virtual_env. Poofs the changes to the virtual environment without any extra cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
