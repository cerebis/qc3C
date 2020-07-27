FROM jcrist/alpine-conda:latest as builder

ENV PYTHONDONTWRITEBYTECODE=true

ENV PATH=/opt/conda/bin:/bin:/usr/bin:/usr/local/bin

RUN conda install --yes --freeze-installed \
	  -c cerebis -c conda-forge -c bioconda \
	  jellyfish \
	  leven \
	  git \
	  && conda clean -afy \
	  && pip --no-cache-dir install git+https://github.com/cerebis/qc3C \
	  && find /opt/conda/ -follow -type f -name '*.a' -delete \
      && find /opt/conda/ -follow -type f -name '*.pyc' -delete

CMD ["/bin/ash"]
