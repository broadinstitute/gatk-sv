FROM alpine:latest
RUN apk --no-cache add curl && \
    curl -s https://api.github.com/repos/Illumina/ExpansionHunterDenovo/releases/latest \
    | grep "browser_download_url.*.tar.gz" \
    | cut -d : -f 2,3 \
    | tr -d \" \
    | xargs wget && \
    mkdir ehdn_extract && \
    tar -xf *.tar.gz --strip-components=1 -C ehdn_extract && \
    rm -rf *.tar.gz && \
    mkdir ehdn && \
    mv ehdn_extract/bin/* ehdn/ && \
    rm -rf ehdn_extract
ENV PATH="/ehdn/:$PATH"
