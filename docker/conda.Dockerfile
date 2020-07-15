FROM continuumio/miniconda3:latest

RUN conda install --yes \
	-c cerebis -c conda-forge \
	 jellyfish \
	 tini

ENTRYPOINT ["tini", "-v", "--"]

CMD ["jellyfish", "-h"]
