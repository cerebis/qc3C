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
        bzip2-devel \
        centos-release-scl \
        gcc \
        gcc-c++ \
        gettext \
        git \
        libtool \
        make \
        ncurses-devel \
        xz-devel \
        zlib-devel \
        wget && \
    yum install -y rh-python36 devtoolset-6-gcc-c++ && \
    yum clean all -y

# setup more current python and gcc support
SHELL ["/usr/bin/scl", "enable",  "rh-python36", "devtoolset-6"]

# install qc3C
RUN pip install -U pip && \
    LC_CTYPE="en_US.UTF-8" pip3 install git+https://github.com/cerebis/qc3C

# jellyfish
WORKDIR /usr/local
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz &&  tar xzf jellyfish-2.2.10.tar.gz
WORKDIR /usr/local/jellyfish-2.2.10
RUN ./configure --prefix=/usr/local && make -j4 && make install
WORKDIR /usr/local/jellyfish-2.2.10/swig/python
RUN export PKG_CONFIG_PATH=/usr/local/jellyfish-2.2.10 && python3 setup.py build && python3 setup.py install
# htslib
WORKDIR /usr/local
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && tar xjf htslib-1.9.tar.bz2
WORKDIR /usr/local/htslib-1.9
RUN make -j4
# samtools
WORKDIR /usr/local/
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xjf samtools-1.9.tar.bz2
WORKDIR /usr/local/samtools-1.9
RUN sed -i -r 's/^(LDFLAGS +=)"?(.*)"?/\1\2 -ltinfo/' Makefile && make -j4
# bwa
WORKDIR /usr/local
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && tar xjf bwa-0.7.17.tar.bz2
WORKDIR bwa-0.7.17
RUN make
# spades
WORKDIR /usr/local
RUN wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz && tar xzf SPAdes-3.13.0-Linux.tar.gz
# bbmap
WORKDIR /usr/local
RUN wget -O bbmap.tar.gz "https://downloads.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2Ffiles%2Flatest%2Fdownload&ts=1555391525" && \
    tar xzf bbmap.tar.gz

#
# stage two - runtime
#
FROM centos:centos6

# create up to date system
RUN yum update -y && \
    yum install -y https://dl.fedoraproject.org/pub/epel/epel-release-latest-6.noarch.rpm && \
    yum install -y --setopt=tsflags=nodocs \
        bzip2 \
        centos-release-scl \
        gettext \
        ncurses \
        pigz \
        xz \
        zlib && \
    yum install -y rh-python36 && \
    yum clean all -y && \
    ldconfig

# copy over the installed software
WORKDIR /usr/local
COPY --from=builder /opt/rh/rh-python36/root /opt/rh/rh-python36/root/
COPY --from=builder /usr/local/bin/jellyfish bin/
COPY --from=builder /usr/local/lib/libjellyfish* lib/
COPY --from=builder /usr/local/lib/pkgconfig lib/pkgconfig/
COPY --from=builder /usr/local/samtools-1.9/samtools /usr/local/bin/
COPY --from=builder /usr/local/htslib-1.9/bgzip /usr/local/bin/
COPY --from=builder /usr/local/bwa-0.7.17/bwa /usr/local/bin/
COPY --from=builder /usr/local/bbmap /usr/local/bbmap
COPY --from=builder /usr/local/SPAdes-3.13.0-Linux /usr/local/SPAdes-3.13.0-Linux/
RUN ln -s /usr/local/bbmap/*sh /usr/local/bin && \
    ln -s /usr/local/SPAdes-3.13.0-Linux/bin/* /usr/local/bin/ && \
    chmod 755 /opt/rh/rh-python36/root/usr/lib/python3.6/site-packages

# further setup
COPY ./root/ /

ENV HOME=/opt/app-root/ \
    PATH=/opt/app-root/bin:/opt/rh/rh-python36/root/usr/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mkdir -p ${HOME}

WORKDIR ${HOME}

ENV BASH_ENV=/opt/app-root/etc/scl_enable \
    ENV=/opt/app-root/etc/scl_enable \
    PROMPT_COMMAND=". /opt/app-root/etc/scl_enable"

# Set the default CMD to print the usage of the language image
ENTRYPOINT ["container-entrypoint"]
CMD ["usage"]
