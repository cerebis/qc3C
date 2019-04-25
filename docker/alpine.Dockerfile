FROM python:3-alpine as builder
MAINTAINER Matthew DeMaere "matt.demaere@gmail.com"

# create up to date system
RUN apk add --no-cache --update-cache \
  automake \
  autoconf \
  build-base \
  bzip2-dev \
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

# install jellyfish binaries and python hooks
WORKDIR /usr/local
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz && tar xzf jellyfish-2.2.10.tar.gz
WORKDIR /usr/local/jellyfish-2.2.10
RUN ./configure && make && make install
WORKDIR /usr/local/jellyfish-2.2.10/swig/python
RUN export PKG_CONFIG_PATH=/usr/local/jellyfish-2.2.10 && python3 setup.py build && python3 setup.py install
WORKDIR /home

FROM python:3-alpine
WORKDIR /usr/local
COPY --from=builder /bin/bash /bin/
COPY --from=builder /usr/lib /usr/lib/
COPY --from=builder /usr/local/bin/jellyfish /usr/local/bin/qc3C /usr/local/bin/tqdm bin/
COPY --from=builder /usr/local/lib/libjellyfish* lib/
COPY --from=builder /usr/local/lib/pkgconfig lib/pkgconfig/
COPY --from=builder /usr/local/lib/python3.7/site-packages lib/python3.7/site-packages/

# further setup
COPY ./root/ /

ENV HOME=/opt/app-root/ \
    PATH=/opt/app-root/bin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p ${HOME}

WORKDIR ${HOME}

# Set the default CMD to print the usage of the language image
ENTRYPOINT ["/usr/bin/container-entrypoint"]
CMD ["usage"]