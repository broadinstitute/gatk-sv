# This image can run codes from the svtk and sv-pipeline repos

# IMPORTANT: these arguments must be specified at the beginning to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_IMAGE=gatksv/sv-pipeline:latest
FROM ${SV_PIPELINE_IMAGE}

RUN pip install pysam==0.15.4 && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends wget make cmake automake && \
    wget -q https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xjf htslib-1.13.tar.bz2 && \
    cd htslib-1.13 && \
    ./configure --prefix=/opt/htslib && \
    make && \
    make install && \
    cd .. && \
    rm -r htslib-1.13* && \
    wget -q https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar -xjf bcftools-1.13.tar.bz2 && \
    cd bcftools-1.13 && \
    ./configure --prefix=/opt/bcftools && \
    make && \
    make install && \
    cd .. && \
    rm -r bcftools-1.13* && \
    apt-get -qqy purge wget make cmake automake && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
ENV PATH="/opt/htslib/bin:/opt/bcftools/bin:${PATH}"
