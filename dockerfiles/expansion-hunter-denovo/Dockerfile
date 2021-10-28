FROM python:3.7-slim
RUN apt-get update && apt-get install -y \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* && \

    wget https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz && \
    mkdir ehdn_extract && \
    tar -xf *.tar.gz --strip-components=1 -C ehdn_extract && \
    rm -rf *.tar.gz && \
    mkdir ehdn && \
    mv ehdn_extract/bin/* ehdn/ && \
    mv ehdn_extract/scripts ehdn/ && \
    rm -rf ehdn_extract && \
    pip install -r /ehdn/scripts/requirements.txt
ENV PATH="/ehdn/:$PATH"
ENV SCRIPTS_DIR /ehdn/scripts
