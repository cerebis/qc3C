FROM python:3-alpine as builder
MAINTAINER Matthew DeMaere "matt.demaere@gmail.com"

# create up to date system
RUN apk add --no-cache --update-cache \
  automake \
  autoconf \
  build-base \
  bzip2-dev \
  ncurses-dev \
  gcc \
  gettext \
  gfortran \
  git \
  libtool \
  openblas-dev \
  wget \
  xz-dev \
  zlib-dev

# install qc3C and dependencies in required order
RUN pip3 install cython && \
  pip3 install numpy && \
  pip3 install git+https://github.com/cerebis/qc3C

# install supporting binaries and python hooks
WORKDIR /usr/local
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz && tar xzf jellyfish-2.2.10.tar.gz
WORKDIR /usr/local/jellyfish-2.2.10
RUN ./configure && make && make install
WORKDIR /usr/local/jellyfish-2.2.10/swig/python
RUN export PKG_CONFIG_PATH=/usr/local/jellyfish-2.2.10 && python3 setup.py build && python3 setup.py install
WORKDIR /home
WORKDIR /usr/local
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar xjf htslib-1.9.tar.bz2
WORKDIR /usr/local/htslib-1.9
RUN make
WORKDIR /usr/local/
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xjf samtools-1.9.tar.bz2
WORKDIR /usr/local/samtools-1.9
RUN make

FROM python:3-alpine
WORKDIR /usr/local
COPY --from=builder /bin/bash /bin/
COPY --from=builder /usr/lib /usr/lib/
COPY --from=builder /usr/local/bin/jellyfish /usr/local/bin/qc3C /usr/local/bin/tqdm bin/
COPY --from=builder /usr/local/lib/libjellyfish* lib/
COPY --from=builder /usr/local/lib/pkgconfig lib/pkgconfig/
COPY --from=builder /usr/local/lib/python3.7/site-packages lib/python3.7/site-packages/
COPY --from=builder /usr/local/samtools-1.9/samtools /usr/local/bin/
COPY --from=builder /usr/local/htslib-1.9/bgzip /usr/local/bin/

# further setup
COPY ./root/usr /usr/

ENV HOME=/opt/app-root/ \
    PATH=/opt/app-root/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p ${HOME}

WORKDIR ${HOME}

# Set the default CMD to print the usage of the language image
ENTRYPOINT ["/usr/bin/container-entrypoint"]
CMD ["usage"]
