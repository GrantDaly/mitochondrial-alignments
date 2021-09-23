FROM debian:bullseye

ENV PATH /root/miniconda3/bin:$PATH
RUN apt-get update -y && apt-get install -y wget && \
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh && \
bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda config --set auto_activate_base true && \
conda install -y bwa && \
conda install -y cutadapt && apt-get remove -y wget


RUN apt-get update -y  && \
apt-get upgrade -y && \
apt-get install -y libgsl-dev pkg-config \
autoconf automake make \
perl zlib1g-dev libbz2-dev liblzma-dev \
libcurl4-gnutls-dev libssl-dev libncurses5-dev \
libperl-dev libgsl0-dev wget gcc g++ git && \
cd ~ && git clone --recursive https://github.com/samtools/htslib.git && \
cd htslib && autoreconf -i && sh configure && make && make install && \
cd ~ &&  git clone --recursive https://github.com/samtools/samtools.git && \
cd samtools && autoreconf && sh configure && make && make install && \
cd ~ && git clone --recursive https://github.com/GrantDaly/bcftools.git && \
cd bcftools && autoheader && autoconf && ./configure --enable-libgsl && \
make && make install && \
cd ~ && rm -r htslib samtools bcftools

ADD cpp /cpp/
RUN ldconfig && \
g++ --std=c++20 -g -Wall -O3 -o inserts-and-cov /cpp/inserts-and-cov.cpp \
`pkg-config --cflags gsl htslib` `pkg-config --libs gsl htslib` && rm -r cpp/ && \
cp inserts-and-cov /usr/local/bin

