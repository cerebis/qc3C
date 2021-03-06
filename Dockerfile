FROM continuumio/miniconda3:latest

LABEL maintainer="matt.demaere@gmail.com"
LABEL org.label-schema.schema-version="1.0"
LABEL org.label-schema.name="cerebis/qc3c"
LABEL org.label-schema.description="qc3C - reference-free quality control for Hi-C DNA sequencing libraries"
LABEL org.label-schema.url="http://github.com/cerebis/qc3C/"
LABEL org.label-schema.vcs-url="http://github.com/cerebis/qc3C/"
LABEL org.label-schema.vcs-ref="e8bca3c1ffbfd15a48e8c784802efc1459719725"
LABEL org.label-schema.version="0.5rc6"
LABEL org.label-schema.docker.cmd="docker run -v /path/to/data:/app cerebis/qc3c kmer -y -m 210 -e DpnII -r reads.fq.gz"

RUN conda install --yes -c cerebis -c conda-forge -c bioconda qc3c "python<3.8" && conda clean -afy

RUN mkdir -p /app
WORKDIR /app
ENTRYPOINT ["qc3C"]
CMD ["--help"]
