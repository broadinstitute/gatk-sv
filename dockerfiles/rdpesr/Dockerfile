# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

# Start with linux base
FROM ubuntu:18.04

# Set environment variables.
ENV HOME=/root
ENV R_HOME=/usr

# Define default command.
CMD ["bash"]

#### Basic image utilities

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get update && apt-get install -y r-base
RUN apt-get update --fix-missing && \
    apt-get upgrade -y && \
    apt-get install -y \
    bcftools \
    bedtools \
    openjdk-8-jdk \
    xvfb \
    zip \
    python3-pip \
    python3-dev \
    curl \
    git \
    samtools \
    tabix \
    vcftools \
    wget \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get -y clean && \
    apt-get -y autoremove && \
    apt-get  -y autoclean

RUN pip3 install --no-cache-dir numpy && \
    pip3 install pysam



# R packages
ARG DEBIAN_FRONTEND=noninteractive
ARG MNORMT_URL="https://cran.r-project.org/src/contrib/Archive/mnormt/mnormt_1.5-7.tar.gz"
ARG MNORMT_DEST="/tmp/mnormt.tar.gz"
ARG DAAG_URL="https://cran.r-project.org/src/contrib/Archive/DAAG/DAAG_1.22.tar.gz"
ARG DAAG_DEST="/tmp/DAAG.tar.gz"
ARG SV_PIPELINE_BIOC_PKGS="\"multtest\""
ARG SV_PIPELINE_R_PKGS="BSDA caret data.table fpc hash metap perm plyr pwr reshape ROCR zoo"
ARG SLIM_R_LIB_CMD="find . -type d \\( -name \"help\" -o -name \"doc\" -o -name \"html\" -o -name \"htmlwidgets\" -o -name \"demo\" -o -name \"demodata\" -o -name \"examples\" -o -name \"exampleData\" -o -name \"unitTests\" -o -name \"tests\" -o -name \"testdata\" -o -name \"shiny\" \\) | xargs rm -rf"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake \
                 libssh2-1-dev \
                 libssl-dev && \
    Rscript -e "install.packages(c('BiocManager', 'latticeExtra','optparse'), repos = 'https://cran.rstudio.com', quiet = TRUE)" && \
    Rscript -e "BiocManager::install(${SV_PIPELINE_BIOC_PKGS})" && \
    curl "${MNORMT_URL}" --output "${MNORMT_DEST}" && \
    curl "${DAAG_URL}" --output "${DAAG_DEST}" && \
    Rscript -e "install.packages(c('${MNORMT_DEST}','${DAAG_DEST}'), repos = NULL, quiet = TRUE)" && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy purge make cmake automake && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

VOLUME ["/root/.config"]

COPY modify_bed_for_PE_SR_RD_labeling.R /src/modify_bed_for_PE_SR_RD_labeling.R
COPY add_SR_PE_to_PB_INS.V2.py /src/add_SR_PE_to_PB_INS.V2.py
COPY add_RD_to_SVs.py     /src/add_RD_to_SVs.py
COPY bincov_to_normCov.R  /src/bincov_to_normCov.R
COPY add_GC_anno_to_bed.R /src/add_GC_anno_to_bed.R
COPY Modify_vcf_by_steps.py /src/Modify_vcf_by_steps.py
COPY integrate_annotations.R /src/integrate_annotations.R


