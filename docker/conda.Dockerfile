FROM continuumio/miniconda3:latest

RUN conda install --yes \
	-c cerebis -c conda-forge -c bioconda \
	 qc3c \
	 samtools \
	 bwa \
	 fastp \
	 tini \
	 && conda clean -afy

ENTRYPOINT ["tini", "--"]

CMD ["qc3C", "-h"]
