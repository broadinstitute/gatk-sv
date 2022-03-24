ARG SVBASE_IMAGE=sv-base:latest
ARG VIRTUAL_ENV_IMAGE=cnmops-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN  rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

# Start from base with most R/Python dependencies included
FROM $SVBASE_IMAGE

# Copy in updated virtual environment
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH

# Add WGD
COPY src/WGD /opt/WGD
ENV PATH="/opt/WGD/bin:${PATH}"

# install supporting libraries
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 file g++ gfortran liblapack-dev libopenblas-dev libxml2-dev && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
