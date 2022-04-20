ARG VIRTUAL_ENV_IMAGE=sv-pipeline-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as build_virtual_env

ARG DEBIAN_FRONTEND=noninteractive
ARG SV_PIPELINE_R_PKGS="assertthat BH cli crayon digest dplyr fansi generics glue magrittr MASS Matrix nlme nnet numDeriv pillar pkgconfig plogr plyr purrr R6 RColorBrewer Rcpp reshape2 rlang rpart stringi stringr survival tibble tidyr tidyselect utf8 withr zoo"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 make cmake automake \
                 libssh2-1-dev \
                 libssl-dev && \
    install_R_packages.R ${SV_PIPELINE_R_PKGS}


# Copy packages from build_virtual_env. Poofs the changes to the virtual environment without any extra cruft.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM $VIRTUAL_ENV_IMAGE
COPY --from=build_virtual_env $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
