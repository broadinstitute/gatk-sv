# This image can run codes from the svtk and sv-pipeline repos

# IMPORTANT: these arguments must be specified at the beginning to take advantage of multi-stage build AND runtime specification of base images
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-virtual-env:latest

FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

# update pysam
RUN pip install pysam==0.15.4

# get build tools
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
        wget make cmake automake \
        libcurl4-openssl-dev libssl-dev

# build htlsib 1.13
RUN wget -q https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure --quiet --prefix=$VIRTUAL_ENV_PATH/htslib && \
    make -s && \
    make -s install

# build bcftools 1.13
RUN wget -q https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar -xjf bcftools-1.13.tar.bz2 && \
    cd bcftools-1.13 && \
    ./configure --quiet --prefix=$VIRTUAL_ENV_PATH/bcftools && \
    make -s && \
    make -s install

# Copy python packages from build_virtual_env. Ppoofs the changes to the virtual environment without any apt cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
ENV PATH=$VIRTUAL_ENV_PATH/htslib/bin:$VIRTUAL_ENV_PATH/bcftools/bin:$PATH
