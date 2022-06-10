ARG SVBASE_IMAGE=sv-base:latest
ARG VIRTUAL_ENV_IMAGE=cnmops-virtual-env:latest

################# get alias for updating dependencies
FROM $SVBASE_IMAGE as base

################ remove unneeded R library files and update required apt dependencies before copying
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh

RUN mv $APT_REQUIRED_PACKAGES /tmp/venv-packages.list
copy --from=base $APT_REQUIRED_PACKAGES $APT_REQUIRED_PACKAGES

RUN export NEW_PACKAGES=$(diff_of_lists.sh /tmp/venv-packages.list $APT_REQUIRED_PACKAGES) && \
    printf " $NEW_PACKAGES" | fix_spaces.sh >> $APT_REQUIRED_PACKAGES

################# Start from base with most R/Python dependencies included
FROM $SVBASE_IMAGE

# Copy in updated virtual environment
COPY --from=virtual_env_image /opt /opt

# Add WGD
COPY src/WGD /opt/WGD
ENV PATH="/opt/WGD/bin:$PATH"

# install any missing support libraries
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends $(fix_spaces.sh $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

RUN runcnMOPS.R  --help