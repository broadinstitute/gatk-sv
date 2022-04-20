# GATK SV Pipeline RdTest dockerfile
ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-rdtest-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

FROM $SV_PIPELINE_IMAGE

# Install needed libraries
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 fftw3 fftw3-dev libssh2-1-dev libssl-dev libjpeg-dev libxml2-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Copy updated virtual environment from VIRTUAL_ENV_IMAGE
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH
