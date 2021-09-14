FROM continuumio/miniconda3

RUN conda config --add channels defaults && \
 conda config --add channels bioconda && \
conda config --add channels conda-forge && \
 conda install bwa && conda update bwa && \
conda install samtools && conda update samtools && \
conda install cutadapt && conda update cutadapt
