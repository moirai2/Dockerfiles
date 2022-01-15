FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
bwa \
cutadapt \
samtools \
bedtools \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
ADD bin /usr/local/bin
ADD config /etc/heliscope
