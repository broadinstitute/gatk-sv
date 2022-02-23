# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_IMAGE=gatksv/sv-pipeline-base:latest
FROM ${SV_PIPELINE_BASE_IMAGE}

# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG SV_PIPELINE_R_PKGS="assertthat BH cli crayon digest dplyr fansi generics glue magrittr MASS Matrix nlme nnet numDeriv pillar pkgconfig plogr plyr purrr R6 RColorBrewer Rcpp reshape2 rlang rpart stringi stringr survival tibble tidyr tidyselect utf8 withr zoo"
ARG SLIM_R_LIB_CMD="find . -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake \
                 libssh2-1-dev \
                 libssl-dev && \
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
