# GATK SV Pipeline RdTest dockerfile

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_R_IMAGE=gatksv/sv-pipeline-children-r:latest
FROM ${SV_PIPELINE_BASE_R_IMAGE}

# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG MNORMT_URL="https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
ARG MNORMT_DEST="/tmp/mnormt.tar.gz"
ARG DAAG_URL="https://cran.r-project.org/src/contrib/Archive/DAAG/DAAG_1.22.tar.gz"
ARG DAAG_DEST="/tmp/DAAG.tar.gz"
ARG SV_PIPELINE_BIOC_PKGS="\"multtest\""
ARG SV_PIPELINE_R_PKGS="BSDA caret fpc hash metap perm plyr pwr reshape ROCR zoo"
ARG SLIM_R_LIB_CMD="find . -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake \
                 libssh2-1-dev \
                 libssl-dev && \
    /opt/install_deprecated_R_package.sh "https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.12.8.tar.gz" && \
    Rscript -e "install.packages(c('BiocManager', 'latticeExtra'), repos = 'https://cran.rstudio.com', quiet = TRUE)" && \
    Rscript -e "BiocManager::install(${SV_PIPELINE_BIOC_PKGS})" && \
    curl "${MNORMT_URL}" --output "${MNORMT_DEST}" && \
    curl "${DAAG_URL}" --output "${DAAG_DEST}" && \
    Rscript -e "install.packages(c('${MNORMT_DEST}','${DAAG_DEST}'), repos = NULL, quiet = TRUE)" && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
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
