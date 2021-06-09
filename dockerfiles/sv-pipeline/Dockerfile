# This image can run codes from the svtk and sv-pipeline repos

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_IMAGE=gatksv/sv-pipeline-base:latest
FROM ${SV_PIPELINE_BASE_IMAGE}

# Install R & bioconductor packages
ARG SV_PIPELINE_BIOCONDUCTOR_PKGS="SNPRelate"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
    Rscript --vanilla install_bioconductor_packages.R ${SV_PIPELINE_BIOCONDUCTOR_PKGS} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy remove --purge make cmake automake && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    bcftools --help

# Install plink2 & KING (for relatedness inference)
ARG KING_URL="https://www.kingrelatedness.com/executables/Linux-king215.tar.gz"
ARG PLING2_URL="https://github.com/chrchang/plink-ng/releases/download/2019/plink2_linux_x86_64_20190107.zip"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 wget unzip && \
    cd /opt && \
    wget -q ${KING_URL} && \
    tar -xzf Linux-king215.tar.gz && \
    wget -q ${PLING2_URL} && \
    unzip plink2_linux_x86_64_20190107.zip && \
    apt-get -qqy purge wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
ENV PATH="/opt/:${PATH}"
