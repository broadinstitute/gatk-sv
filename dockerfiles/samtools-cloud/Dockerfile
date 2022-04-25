# This docker file exists because we need samtools to be able to access GCS buckets
# so we base off the sv-base-mini image, which has samtools installed, 
# and only install GCloud SDK here for authentication purpose

# Start with the barebones image that has samtools (a version that must support NIO) installed 
ARG MINIBASE_IMAGE=sv-base-mini:latest
ARG VIRTUAL_ENV_IMAGE=samtoolscloud-virtual-env:latest
# available gsutil versions here: https://cloud.google.com/sdk/docs/release-notes
ARG CLOUD_SDK_VERSION=389.0.0
ARG CONDA_ENV_NAME="gatk-sv"
ARG CONDA_INSTALL_DIR=/opt/conda
ARG GOOGLE_CLOUD_INSTALL_DIR=/opt/google-cloud-sdk

################## install google-cloud-cli into a copy of the virtual image
FROM $VIRTUAL_ENV_IMAGE as prep_virtual_env

# get prerequisites for installing google-cloud-cli
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 apt-transport-https ca-certificates wget libssl3

# install google-cloud-cli
# NOTE --override-components to avoid installing bq (BigQuery)
# NOTE: despite this being a "versioned" package, the installer will insist on upgrading at install time, so this
#       is difficult to fully pin
ARG CLOUD_SDK_VERSION
ARG GOOGLE_CLOUD_INSTALL_DIR
ENV CLOUDSDK_CONFIG=$GOOGLE_CLOUD_INSTALL_DIR
ENV PATH=$GOOGLE_CLOUD_INSTALL_DIR/bin:$PATH
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${CLOUD_SDK_VERSION}-linux-x86_64.tar.gz && \
    tar xzf google-cloud-sdk-${CLOUD_SDK_VERSION}-linux-x86_64.tar.gz && \
    mkdir -p $(basename $GOOGLE_CLOUD_INSTALL_DIR) && \
    mv google-cloud-sdk $GOOGLE_CLOUD_INSTALL_DIR && \
    cd $GOOGLE_CLOUD_INSTALL_DIR && \
    ./install.sh --usage-reporting false --override-components core gsutil && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image

# clean up some stuff we don't need
#    -anthoscli  (Anthos)
RUN rm $GOOGLE_CLOUD_INSTALL_DIR/bin/anthoscli

########## Copy results to final image
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

# show the google packages run
RUN gcloud --version
RUN gsutil --version
