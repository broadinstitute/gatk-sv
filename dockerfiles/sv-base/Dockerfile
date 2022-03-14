# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, RdTest, and VCF-QC images
#

# Start with the barebones image
ARG MINIBASE_IMAGE=gatksv/sv-base-mini:latest
FROM ${MINIBASE_IMAGE}

ENV HOME=/root

####################################################################################################
# GATK needed for PrintSVEvidence
ARG GATK_TRANSIENT_DEP="git git-lfs"
ARG GATK_DEP="openjdk-8-jdk ${GATK_TRANSIENT_DEP}"
ARG GATK_COMMIT="a33bf19dd3188af0af1bd17bce015eb20ba73227"
ENV GATK_JAR="/opt/gatk.jar"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${GATK_DEP} && \
    cd / && git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && git lfs install && \
    git checkout ${GATK_COMMIT} && ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR} && \
    cd / && \
    apt-get -qqy purge ${GATK_TRANSIENT_DEP} && \
    apt-get -qqy clean && \
    rm -rf /gatk \
           /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    java -jar ${GATK_JAR}

# we settle on this order because the autoremove steps in G-Cloud and R blocks simply hate python and remove it

####################################################################################################
# G-Cloud
ARG CLOUD_SDK_VERSION=235.0.0
ARG GCLOUD_SDK_EXTRA_COMPONENTS
ARG GCLOUD_SDK_TRANSIENT_DEP="gpg-agent gnupg python-pip python-setuptools "
ARG GCLOUD_SDK_DEPENDENCIES="lsb-release ${GCLOUD_SDK_TRANSIENT_DEP}"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${GCLOUD_SDK_DEPENDENCIES} && \
    pip install -q -U --no-cache-dir \
                crcmod \
                wheel && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 google-cloud-sdk=${CLOUD_SDK_VERSION}-0 \
                 $GCLOUD_SDK_EXTRA_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy purge ${GCLOUD_SDK_TRANSIENT_DEP} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    bcftools --help

####################################################################################################
# R and all packages
# Steps:
#  1. OS libs
#  2. install base R, see https://cran.r-project.org/bin/linux/ubuntu/README.html#administration-and-maintances-of-r-packages
#  3. install R packages that all/most child images depend on
#  4. remove doc, help, demos and tests (saves ~ 10% space)
#  4. final clean up
COPY dockerfiles/sv-base/install_R_packages.R /opt/
COPY dockerfiles/sv-base/install_bioconductor_packages.R /opt/
COPY dockerfiles/sv-base/install_deprecated_R_package.sh /opt/
ARG R_DEP_TRANSIENT="make gpg-agent gnupg"
ARG R_DEPENDENCIES="software-properties-common build-essential \
    ${R_DEP_TRANSIENT} libz-dev libncurses5-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev ca-certificates \
    apt-transport-https openssh-client"
ARG R_RELEASE_VERSION="3.5.1"
ARG SV_BASE_R_PKGS="optparse BiocManager"
ARG SLIM_R_LIB_CMD="find .  -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${R_DEPENDENCIES} && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/" && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 r-base-core="${R_RELEASE_VERSION}"-2bionic \
                 r-base-dev="${R_RELEASE_VERSION}"-2bionic && \
    mkdir -p /tmp/R_pkg_download/ && \
    Rscript -e "update.packages(destdir = \"/tmp/R_pkg_download/\", clean = TRUE, quiet = TRUE)" && \
    cd /opt/ && \
    Rscript --vanilla install_R_packages.R ${SV_BASE_R_PKGS} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy purge \
                 ${R_DEP_TRANSIENT} && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
