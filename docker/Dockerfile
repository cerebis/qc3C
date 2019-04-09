FROM fedora:29
MAINTAINER Matthew DeMaere "matt.demaere@gmail.com"

# create up to date system
RUN dnf update -y && dnf install -y \
        autoconf \
        automake \
        gcc \
        gcc-c++ \
        gettext \
        git \
        libtool \
        make \
        python3 \
        python3-devel \
        python3-pip \
        wget

# install qc3C
RUN pip3 install -U pip && pip3 install git+https://github.com/cerebis/qc3C

# install jellyfish binaries and python hooks
WORKDIR "/usr/local"
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz && tar xzf jellyfish-2.2.10.tar.gz
WORKDIR "/usr/local/jellyfish-2.2.10"
RUN ./configure && make && make install
WORKDIR "/usr/local/jellyfish-2.2.10/swig/python"
RUN export PKG_CONFIG_PATH=/usr/local/jellyfish-2.2.10 && python3 setup.py build && python3 setup.py install
WORKDIR /home
