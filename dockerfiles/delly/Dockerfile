# use the ubuntu base image
FROM ubuntu:14.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    wget \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr
ENV EBROOTHTSLIB /opt/htslib

# install delly, pinned to a particular commit
ARG HTSLIB_VERSION=1.9
RUN cd /opt \
    && wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && rm htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && git checkout 65863d9aef37a897b1a0116bd7c992ab64755f88 \
    && make STATIC=1 all \
    && make install

# Multi-stage build
FROM alpine:latest
RUN apk add --no-cache bash
RUN mkdir -p /opt/delly/bin
WORKDIR /opt/delly/bin
COPY --from=0 /opt/delly/bin/delly .

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

WORKDIR /

# by default /bin/sh is executed
CMD ["/bin/bash"]
