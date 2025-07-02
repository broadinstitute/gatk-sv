# This docker file exists because we need samtools to be able to access GCS buckets
# so we base off the sv-base-mini image, which has samtools installed, 
# and only install GCloud SDK here for authentication purposes.

# Start with the barebones image that has samtools (a version that must support NIO) installed
ARG MINIBASE_IMAGE=sv-base-mini:latest
ARG VIRTUAL_ENV_IMAGE=samtoolscloud-virtual-env:latest
# available gsutil versions here: https://cloud.google.com/sdk/docs/release-notes
ARG CLOUD_SDK_VERSION=518.0.0-0
ARG CONDA_ENV_NAME="gatk-sv"
ARG CONDA_INSTALL_DIR=/opt/conda
################## install google-cloud-cli into a copy of the virtual image
FROM $VIRTUAL_ENV_IMAGE as prep_virtual_env
# get prerequisites for installing google-cloud-cli
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 apt-transport-https ca-certificates wget curl gnupg2 libssl3 \
                 libfontconfig1-dev

# Google cloud SDK is installed using the documentation in the following link.
# https://cloud.google.com/sdk/docs/install#deb
ARG CLOUD_SDK_VERSION

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && \
    apt-get install google-cloud-cli=$CLOUD_SDK_VERSION -y && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image

ENV GOOGLE_CLOUD_INSTALL_DIR=/usr/lib/google-cloud-sdk
ENV CLOUDSDK_CONFIG=$GOOGLE_CLOUD_INSTALL_DIR

# clean up some stuff we don't need
#    -anthoscli  (Anthos)
RUN rm $GOOGLE_CLOUD_INSTALL_DIR/bin/anthoscli
# Copy results to final image
FROM $MINIBASE_IMAGE
ARG CONDA_INSTALL_DIR
ARG CONDA_ENV_NAME
ARG CONDA_ENV_PATH=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME
# copy conda environment
COPY --from=prep_virtual_env $CONDA_ENV_PATH $CONDA_ENV_PATH
# copy google-cloud-cli
ARG GOOGLE_CLOUD_INSTALL_DIR
COPY --from=prep_virtual_env $GOOGLE_CLOUD_INSTALL_DIR $GOOGLE_CLOUD_INSTALL_DIR
ENV PATH=/opt/google-cloud-sdk/bin:$CONDA_ENV_PATH/bin:$PATH
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends libfontconfig1-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
# show the google packages run
RUN gcloud --version
RUN gsutil --version
