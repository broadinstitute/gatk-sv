# GATK SV Pipeline Hail dockerfile

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_IMAGE=gatksv/sv-pipeline:latest
FROM ${SV_PIPELINE_IMAGE}

# Dependencies for creating a Hail cluster on GCS DataProc
RUN pip3 --no-cache-dir install hail==0.2.71 && \
    pip3 --no-cache-dir install google-cloud-dataproc \
