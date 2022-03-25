# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, RdTest, and VCF-QC images
ARG MINIBASE_IMAGE=sv-base-mini:latest
ARG VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG VIRTUAL_ENV_PATH="/gatk-sv-env"

FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
ARG VIRTUAL_ENV_PATH
RUN rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

FROM $MINIBASE_IMAGE
# Start with the barebones image

ENV HOME=/root

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
           /usr/share/man/??_*

# Copy virtual environment and symlink R to its expected location
ARG VIRTUAL_ENV_PATH
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
RUN mkdir -p /opt/R/ &&  \
    export R_RELEASE=$(ls /gatk-sv-env/R/) && \
    ln -s $VIRTUAL_ENV_PATH/R/$R_RELEASE /opt/R/$R_RELEASE && \
    ln -s /opt/R/$R_RELEASE/bin/* $VIRTUAL_ENV_PATH/bin

# Add virtual environment executables to PATH. Now environment is installed and "activated"
ENV PATH=$VIRTUAL_ENV_PATH:$VIRTUAL_ENV_PATH/bin:$PATH
ENV VIRTUAL_ENV_PATH=$VIRTUAL_ENV_PATH
ENV GATK_JAR=$VIRTUAL_ENV_PATH/gatk.jar

# Install libraries to support R and gatk (needed for PrintSVEvidence) and test to ensure it's installed correctly
ARG GATK_DEP="openjdk-8-jdk"
ARG ADDITIONAL_SUPPORT_PACKAGES="file"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    xargs apt-get -qqy install --no-upgrade --no-install-recommends < $VIRTUAL_ENV_PATH/r_support_packages.list && \
    apt-get -qqy install --no-upgrade --no-install-recommends ${GATK_DEP} $ADDITIONAL_SUPPORT_PACKAGES && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    R --version && \
    java -jar ${GATK_JAR} && \
    bcftools --help
