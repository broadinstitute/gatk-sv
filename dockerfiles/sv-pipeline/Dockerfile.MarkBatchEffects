# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This Dockerfile builds the primary image used for most core genomic data 
# processing, filtering, and quality control operations that rely on code from
# this git repo

# It is based on the primary sv-pipeline image from GATK-SV, which is available here:
# us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline
ARG BASE_IMAGE=us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-05-01-v1.0.3-85054e07
FROM $BASE_IMAGE as base

# Aliases for convenience
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc && \
    echo 'alias less="zless"' >> ~/.bashrc

# Install plink2
RUN curl -sSL https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250515.zip -o plink2.zip && \
    unzip plink2.zip -d /usr/local/bin && \
    chmod +x /usr/local/bin/plink2 && \
    rm plink2.zip && \
    plink2 --version

# Install git & lftp
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
        git \
        lftp && \
        apt-get clean

# Install missing R libraries as necessary
RUN Rscript -e 'install.packages("argparse", repos="https://cran.rstudio.com")'

# Clone pancan_germline_wgs repo and check out specified hash
ARG G2C_TAG=main
RUN git clone \
    https://github.com/vanallenlab/pancan_germline_wgs.git \
    opt/pancan_germline_wgs && \
    cd opt/pancan_germline_wgs && \
    git pull && \
    git checkout $G2C_TAG && \
    cd -

# Launch bash as entrypoint
CMD ["/bin/bash"]
