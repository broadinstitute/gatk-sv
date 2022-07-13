ARG SAMTOOLS_CLOUD_IMAGE=samtoolscloud:latest:latest
ARG VIRTUAL_ENV_IMAGE=samtoolscloud-virtual-env:latest

ARG CONDA_ENV_NAME="gatk-sv"
ARG CONDA_INSTALL_DIR=/opt/conda
ARG SV_UTILS_IMAGE_PATH=/opt/sv_utils/
ARG SV_UTILS_LOCAL_PATH=src/sv_utils

# get alias to copy virtual env
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image

########## Copy results to final image
FROM $SAMTOOLS_CLOUD_IMAGE
ARG CONDA_INSTALL_DIR
ARG CONDA_ENV_NAME
ARG CONDA_ENV_PATH=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME
# copy conda environment
COPY --from=virtual_env_image $CONDA_ENV_PATH $CONDA_ENV_PATH
ENV CONDA_ENV_NAME=$CONDA_ENV_NAME
ENV CONDA_INSTALL_DIR=$CONDA_INSTALL_DIR
ARG CONDA_ENV_BIN=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME/bin
ENV PATH=$CONDA_ENV_BIN:$CONDA_BIN:$PATH

# copy sv-utils and installed
COPY src/sv_utils /opt/sv_utils
RUN pip install -e /opt/sv_utils

# show that sv-utils is installed
RUN sv-utils --help
