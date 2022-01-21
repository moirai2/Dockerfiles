#https://cellrank.readthedocs.io/
#https://cellrank.readthedocs.io/en/stable/installation.html#dependencies
FROM python:latest
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV TZ JST-9
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libopenmpi-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
RUN pip install --default-timeout=100  --upgrade --no-cache-dir \
cython \
mpi4py \
petsc \
slepc \
slepc4py \
&& rm -rf /root/.cache/
