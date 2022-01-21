#https://cellrank.readthedocs.io/
#https://cellrank.readthedocs.io/en/stable/installation.html#dependencies
FROM python:latest
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV TZ JST-9
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libopenmpi-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& pip install cython \
&& pip install mpi4py \
&& pip install petsc \
&& pip install slepc \
&& pip install slepc4py \
&& rm -rf /root/.cache/
