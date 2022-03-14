ARG VIRTUAL_ENV_IMAGE=sv-pipeline-base-virtual-env:latest

FROM $VIRTUAL_ENV_IMAGE as build_virtual_env
# get supporting libraries
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake apt-transport-https ca-certificates wget unzip \
                 file g++ gfortran liblapack-dev libopenblas-dev libxml2-dev

# Install R and bioconductor dependencies
ARG SV_PIPELINE_BIOCONDUCTOR_PKGS="SNPRelate"
RUN install_bioconductor_packages.R $SV_PIPELINE_BIOCONDUCTOR_PKGS

# Install plink2 & KING (for relatedness inference)
ARG KING_URL="https://www.kingrelatedness.com/executables/Linux-king215.tar.gz"
RUN cd $VIRTUAL_ENV_PATH && \
    wget -q $KING_URL && \
    tar -xzf Linux-king215.tar.gz && \
    rm -f Linux-king215.tar.gz

ARG PLINK2_URL="https://github.com/chrchang/plink-ng/releases/download/2019/plink2_linux_x86_64_20190107.zip"
RUN cd $VIRTUAL_ENV_PATH && \
    wget -q $PLINK2_URL && \
    unzip plink2_linux_x86_64_20190107.zip && \
    rm -f plink2_linux_x86_64_20190107.zip

# Copy python packages from build_virtual_env. Ppoofs the changes to the virtual environment without any apt cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
