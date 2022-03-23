# This docker file exists because we need samtools to be able to access GCS buckets
# so we base off the sv-base-mini image, which has samtools installed, 
# and only install GCloud SDK here for authentication purpose

# Start with the barebones image that has samtools (a version that must support NIO) installed 
ARG MINIBASE_IMAGE=sv-base-mini:latest
FROM ${MINIBASE_IMAGE}

ENV HOME=/root

# install GCloud SDK
ARG CLOUD_SDK_VERSION=239.0.0
ARG GCLOUD_SDK_EXTRA_COMPONENTS
ARG GCLOUD_SDK_TRANSIENT_DEP="curl gpg-agent gnupg python-pip python-setuptools "
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
