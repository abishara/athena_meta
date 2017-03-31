FROM ubuntu:16.04

RUN apt-get update

RUN apt-get -y install build-essential

RUN apt-get -y install python2.7 && ln -s /usr/bin/python2.7 /usr/bin/python

RUN apt-get -y install openjdk-8-jre

RUN apt-get install wget

RUN apt-get -y install autotools-dev \
     && apt-get -y install automake \
     && apt-get -y install autoconf

RUN apt-get -y install zlib1g-dev

RUN apt-get -y install ncurses-dev

RUN wget https://github.com/grocsvs/idba/archive/1.1.3g1.tar.gz \
     && tar -xf 1.1.3g1.tar.gz \
     && cd idba-1.1.3g1 \
     && ./build.sh \
     && ./configure \
     && make \
     && mv bin/idba_ud /bin

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 \
     && tar -xf bwa-0.7.15.tar.bz2 \
     && cd bwa-0.7.15 \
     && make \
     && mv bwa /bin

RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 \
     && tar -xf samtools-1.3.1.tar.bz2 \
     && cd samtools-1.3.1 \
     && make install

RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 \
     && tar -xf htslib-1.3.2.tar.bz2 \
     && cd htslib-1.3.2 \
     && make install

RUN wget https://github.com/marbl/canu/archive/v1.4.tar.gz \
     && tar -xf v1.4.tar.gz \
     && cd canu-1.4/src \
     && make \
     && ln -s $(readlink -f ../Linux-amd64/bin/canu) /bin/canu

RUN apt-get -y install python-pip && pip install -U pip

RUN mkdir athena_meta_src && cd athena_meta_src \
     && wget https://github.com/abishara/athena_meta/archive/1.0.tar.gz -O athena_meta.tar.gz \
     && tar -xf athena_meta.tar.gz --strip-components 1 \
     && pip install -r requirements.txt \
     && pip install -vvv .

CMD athena-meta
