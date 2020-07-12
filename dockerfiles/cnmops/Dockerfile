# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG GATKSV_PIPELINE_V1_RESOURCES_IMAGE=gatksv-pipeline-v1-resources:latest
ARG SVBASE_IMAGE=gatksv/sv-base:latest

### stage 1
FROM ${GATKSV_PIPELINE_V1_RESOURCES_IMAGE}
RUN mkdir -p /opt/ && \
    mv ${RESOURCES}WGD /opt/WGD

### stage 2
# Start from base with most R/Python dependencies included
FROM ${SVBASE_IMAGE}

# Install R and bioconductor dependencies
ARG SLIM_R_LIB_CMD="find .  -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake \
                 file \
                 g++ \
                 gfortran \
                 liblapack-dev \
                 libopenblas-dev \
                 libxml2-dev && \
    mkdir -p /tmp/R_pkg_download/ && \
    /opt/install_deprecated_R_package.sh https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz && \
    cd /opt/ && \
    Rscript --vanilla install_bioconductor_packages.R cn.mops rtracklayer && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy remove make cmake automake && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Add WGD
COPY --from=0 /opt/WGD /opt/WGD
ENV PATH="/opt/WGD/bin:${PATH}"
