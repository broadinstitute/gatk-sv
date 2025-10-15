# Simple Dockerfile for STRipy-pipeline
# Force x86_64 platform for compatibility with genomics tools
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

# Install minimal dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    wget \
    samtools \
    libcurl4 \
    libssl3 \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install ExpansionHunter (core dependency for STR analysis)
ARG EXPANSIONHUNTER_VERSION=5.0.0
RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/v${EXPANSIONHUNTER_VERSION}/ExpansionHunter-v${EXPANSIONHUNTER_VERSION}-linux_x86_64.tar.gz \
    && tar -xzf ExpansionHunter-v${EXPANSIONHUNTER_VERSION}-linux_x86_64.tar.gz \
    && mv ExpansionHunter-v${EXPANSIONHUNTER_VERSION}-linux_x86_64/bin/ExpansionHunter /usr/local/bin/ \
    && rm -rf ExpansionHunter-v${EXPANSIONHUNTER_VERSION}-linux_x86_64*

# Install REViewer (for read visualization in STR analysis)
ARG REVIEWER_VERSION=0.2.7
RUN wget https://github.com/Illumina/REViewer/releases/download/v${REVIEWER_VERSION}/REViewer-v${REVIEWER_VERSION}-linux_x86_64.gz \
    && gunzip REViewer-v${REVIEWER_VERSION}-linux_x86_64.gz \
    && chmod +x REViewer-v${REVIEWER_VERSION}-linux_x86_64 \
    && mv REViewer-v${REVIEWER_VERSION}-linux_x86_64 /usr/local/bin/REViewer

# Install BWA aligner
ARG BWA_VERSION=0.7.19
RUN wget https://github.com/lh3/bwa/archive/refs/tags/v${BWA_VERSION}.tar.gz \
    && tar -xzf v${BWA_VERSION}.tar.gz \
    && cd bwa-${BWA_VERSION} \
    && make \
    && mv bwa /usr/local/bin/ \
    && cd .. \
    && rm -rf bwa-${BWA_VERSION}* v${BWA_VERSION}.tar.gz

# Clone STRipy-pipeline repository
RUN git clone https://gitlab.com/andreassh/stripy-pipeline.git /opt/stripy-pipeline

# Set working directory
WORKDIR /opt/stripy-pipeline

# Install Python requirements if they exist
RUN if [ -f requirements.txt ]; then pip3 install -r requirements.txt; fi

# Copy wrapper script, default config, and loci profiles
COPY scripts/stripy /usr/local/bin/stripy
COPY scripts/default-config.json /usr/local/bin/default-config.json
COPY scripts/catalogs /opt/stripy-pipeline/scripts/catalogs
COPY scripts/stripy_to_vcf.py /opt/stripy-pipeline/scripts/stripy_to_vcf.py
RUN pip3 install pyfaidx

# Make wrapper executable
RUN chmod +x /usr/local/bin/stripy

# Set entrypoint to use the wrapper
ENTRYPOINT ["/usr/local/bin/stripy"]

# Default help command
CMD ["--help"]
