# use the ubuntu base image
FROM ubuntu:18.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install sansa
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/sansa.git \
    && cd /opt/sansa/ \
    && make STATIC=1 all \
    && make install


# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/sansa/bin
WORKDIR /opt/sansa/bin
COPY --from=0 /opt/sansa/bin/sansa .

# Workdir
WORKDIR /root/

# Add sansa to PATH
ENV PATH="/opt/sansa/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]
