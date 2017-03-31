FROM python:2

RUN apt-get update

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

RUN wget https://github.com/marbl/canu/archive/v1.3.tar.gz \
     && tar -xf v1.3.tar.gz \
     && cd canu-1.3/src \
     && make \
     && ln -s $(readlink -f ../Linux-amd64/bin/canu) /bin/canu

RUN pip install -U pip

CMD bash
