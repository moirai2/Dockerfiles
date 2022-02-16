FROM python:3.7
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
samtools \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
WORKDIR /root
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
pip \
matplotlib \
numba \
numpy \
scipy \
cython \
scikit-learn \
h5py \
click
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
velocyto \
&& rm -rf ~/.cache/pip