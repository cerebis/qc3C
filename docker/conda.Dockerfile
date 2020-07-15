FROM continuumio/miniconda3:latest

RUN conda install --yes \
	-c cerebis -c conda-forge -c bioconda \
	 qc3C \
	 tini

ENTRYPOINT ["tini", "--"]

CMD ["qc3C", "-h"]
