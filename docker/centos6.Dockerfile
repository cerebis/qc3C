FROM centos:centos6 as builder
MAINTAINER Matthew DeMaere "matt.demaere@gmail.com"

ENV SUMMARY="qc3C Hi-C DNA sequencing quality control" \
    DESCRIPTION="qc3C is a tool for assessing the signal strength of Hi-C sequencing read-sets. \
    Assessment can be performed on either the reads alone or against the reads once mapped to an assembly. \
    The results returned are an estimate of the number of Hi-C proximity ligation containing read-pairs."

LABEL name="cerebis/qc3c-centos6" \
      version="0.1" \
      summary="$SUMMARY" \
      description="$DESCRIPTION" \
      usage="docker run -it -v /local/data:/opt/app-root/data:z cerebis/qc3c bash"

# create up to date system
RUN yum update -y && \
    yum install -y --setopt=tsflags=nodocs \
        autoconf \
        automake \
        centos-release-scl \
        gcc \
        gcc-c++ \
        gettext \
        git \
        libtool \
        make \
        wget && \
    yum install -y rh-python36 devtoolset-6-gcc-c++ && \
    yum clean all -y

# setup more current python and gcc support
SHELL ["/usr/bin/scl", "enable",  "rh-python36", "devtoolset-6"]

# install qc3C
RUN pip install -U pip && \
    LC_CTYPE="en_US.UTF-8" pip3 install git+https://github.com/cerebis/qc3C

# install jellyfish
WORKDIR /usr/local

RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz && \
    tar xzf jellyfish-2.2.10.tar.gz

WORKDIR /usr/local/jellyfish-2.2.10

RUN ./configure --prefix=/usr/local && \
    make && \
    make install

# install Python hooks from Jellyfish
WORKDIR /usr/local/jellyfish-2.2.10/swig/python
RUN export PKG_CONFIG_PATH=/usr/local/jellyfish-2.2.10 && \
    python3 setup.py build && \
    python3 setup.py install

# further setup
COPY ./root/ /

ENV HOME=/opt/app-root/ \
    PATH=/opt/app-root/bin:/opt/rh/devtoolset-6/root/usr/bin/:/opt/rh/rh-python36/root/usr/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p ${HOME}

WORKDIR ${HOME}

ENV BASH_ENV=/opt/app-root/etc/scl_enable \
    ENV=/opt/app-root/etc/scl_enable \
    PROMPT_COMMAND=". /opt/app-root/etc/scl_enable"

# Set the default CMD to print the usage of the language image
ENTRYPOINT ["container-entrypoint"]
CMD ["usage"]
