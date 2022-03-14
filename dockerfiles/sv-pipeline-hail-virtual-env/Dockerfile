# GATK SV Pipeline Hail virtual environment dockerfile
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-virtual-env:latest
FROM $VIRTUAL_ENV_IMAGE

# Dependencies for creating a Hail cluster on GCS DataProc
RUN pip3 --no-cache-dir install hail==0.2.71 google-cloud-dataproc
