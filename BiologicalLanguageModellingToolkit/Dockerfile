FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
WORKDIR /root
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
wget \
gcc \
make \
build-essential \
ca-certificates \
&& apt-get autoremove -y \
&& apt-get clean \
&& rm -rf /usr/local/src/*
RUN wget https://www.cs.cmu.edu/~blmt/source/blmt_v1.5.tar.gz \
&& tar -zxf blmt_v1.5.tar.gz \
&& cd blmt_version2 \
&& sed -i 's/CC = CC/CC = gcc/g' makefile \
&& make \
&& mv faa2srt /usr/local/bin/. \
&& mv ngrams /usr/local/bin/. \
&& mv proteinCount /usr/local/bin/. \
&& mv proteinNGram /usr/local/bin/. \
&& mv srt2lcp /usr/local/bin/. \
&& cd .. \
&& rm -r blmt_version2
