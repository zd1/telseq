# telseq docker
#
# To build:
#
# docker build -t telseq-docker github.com/zd1/telseq
#
# To run:
#
# docker run telseq-docker a.bam
#

FROM ubuntu:14.04
MAINTAINER Zhihao Ding <zhihao.ding@gmail.com>
LABEL Description="Telseq docker" Version="0.0.1"

VOLUME /tmp

WORKDIR /tmp

RUN apt-get update && \
    apt-get install -y \
        automake \
        autotools-dev \
        build-essential \
        cmake \
        libhts-dev \
        libhts0 \
        libjemalloc-dev \
        libsparsehash-dev \
        libz-dev \
        python-matplotlib \
        wget \
        zlib1g-dev

# build remaining dependencies:
# bamtools
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xzvf v2.4.0.tar.gz && \
    rm v2.4.0.tar.gz && \
    cd bamtools-2.4.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make


# build telseq
RUN mkdir -p /src && \
    cd /src && \
    wget https://github.com/zd1/telseq/archive/0.0.1.tar.gz && \
    tar -xzvf 0.0.1.tar.gz && \
    rm 0.0.1.tar.gz && \
    cd telseq-0.0.1/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/deps/bamtools-2.4.0 --with-jemalloc=/usr --prefix=/usr/local && \
    make && \
    make install


ENTRYPOINT ["/usr/local/bin/telseq"]
CMD ["--help"]