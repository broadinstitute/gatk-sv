# GATK SV Pipeline VCF-QC / WGD dockerfile

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_R_IMAGE=gatksv/sv-pipeline-children-r:latest
FROM ${SV_PIPELINE_BASE_R_IMAGE}

# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG GERT_DEPENDENCY_R_PKGS="openssl credentials askpass"
ARG SV_PIPELINE_R_PKGS="beeswarm devtools HardyWeinberg nloptr RColorBrewer vioplot zoo"
ARG SLIM_R_LIB_CMD="find . -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake \
                 libgit2-dev \
                 libssh2-1-dev \
                 libxml2-dev \
                 libssl-dev && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
    Rscript --vanilla install_R_packages.R ${GERT_DEPENDENCY_R_PKGS} && \
    /opt/install_deprecated_R_package.sh "https://cran.r-project.org/src/contrib/Archive/gert/gert_0.3.tar.gz" && \
    Rscript --vanilla install_R_packages.R ${SV_PIPELINE_R_PKGS} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy purge make cmake automake && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
