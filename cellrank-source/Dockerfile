#https://cellrank.readthedocs.io/
#https://cellrank.readthedocs.io/en/stable/installation.html#dependencies
#https://www.mcs.anl.gov/petsc/petsc4py-current/docs/usrman/install.html
#https://itneko.blogspot.com/2019/08/ubuntugfortran.html
#https://askubuntu.com/questions/623578/installing-blas-and-lapack-packages
FROM python:latest
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV TZ JST-9
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libopenmpi-dev \
gfortran \
libblas-dev \
liblapack-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
cython \
numpy \
mpi4py \
&& rm -rf /root/.cache/
RUN cd /usr/local \
&& git clone -b release https://gitlab.com/petsc/petsc \
&& cd petsc \
&& ./configure --download-mpich --download-fblaslapack \
&& make all check
RUN cd /usr/local \
&& git clone -b release  https://gitlab.com/slepc/slepc \
&& cd slepc \
&& export PETSC_DIR=/usr/local/petsc \
&& export SLEPC_DIR=/usr/local/slepc \
&& export PETSC_ARCH=arch-linux-c-debug \
&& ./configure \
&& make \
&& make check
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
petsc4py
RUN mkdir -p /usr/local/lib/python3.10/site-packages/petsc/installed-arch-linux2-c-opt/ \
&& ln -s /usr/local/lib/python3.10/site-packages/petsc/include/ /usr/local/lib/python3.10/site-packages/petsc/installed-arch-linux2-c-opt/include \
&& ln -s /usr/local/lib/python3.10/site-packages/petsc/lib/ /usr/local/lib/python3.10/site-packages/petsc/installed-arch-linux2-c-opt/lib
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
slepc4py \
&& rm -rf /root/.cache/
