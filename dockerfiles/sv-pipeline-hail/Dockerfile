FROM us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline-base:rlc_posthoc_filtering_cnv_mcnv_compatability_9a8561
RUN mkdir -p /usr/share/man/man1 && \
    apt-get update && apt-get install -y \
    openjdk-8-jre-headless \
    && rm -rf /var/lib/apt/lists/* && \
    pip3 --no-cache-dir install hail==0.2.71 && \
    pip3 --no-cache-dir install google-cloud-dataproc
ENTRYPOINT ["/bin/bash"]