FROM alpine:latest
RUN apk --no-cache add curl && \
    wget https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz && \
    mkdir ehdn_extract && \
    tar -xf *.tar.gz --strip-components=1 -C ehdn_extract && \
    rm -rf *.tar.gz && \
    mkdir ehdn && \
    mv ehdn_extract/bin/* ehdn/ && \
    rm -rf ehdn_extract
ENV PATH="/ehdn/:$PATH"
