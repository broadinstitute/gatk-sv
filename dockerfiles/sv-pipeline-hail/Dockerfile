# GATK SV Pipeline Hail dockerfile
ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-hail-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

FROM $SV_PIPELINE_IMAGE

# Copy updated virtual environment from VIRTUAL_ENV_IMAGE
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
