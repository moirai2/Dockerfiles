FROM python:3.7
USER root
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV TZ JST-9
RUN apt-get update \
&& apt-get -y install \
vim \
less \
&& apt-get autoremove -y \
&& apt-get clean \
&& rm -rf /usr/local/src/*
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
pip \
jupyterlab \
&& rm -rf ~/.cache/pip
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libopenmpi-dev \
gfortran \
libblas-dev \
liblapack-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
cython==0.29.27 \
numpy==1.21.5 \
mpi4py==3.1.3 \
petsc==3.16.1 \
petsc4py==3.16.1 \
&& rm -rf /root/.cache/
RUN mkdir -p /usr/local/lib/python3.7/site-packages/petsc/installed-arch-linux2-c-opt/ \
&& ln -s /usr/local/lib/python3.7/site-packages/petsc/include /usr/local/lib/python3.7/site-packages/petsc/installed-arch-linux2-c-opt/include \
&& ln -s /usr/local/lib/python3.7/site-packages/petsc/lib /usr/local/lib/python3.7/site-packages/petsc/installed-arch-linux2-c-opt/lib
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
slepc==3.16.1 \
slepc4py==3.16.1 \
&& rm -rf /root/.cache/
COPY requirements.txt /root
WORKDIR /root
RUN pip install --default-timeout=100 --upgrade --no-cache-dir -r requirements.txt \
&& rm requirements.txt \
&& rm -rf ~/.cache/pip
CMD [" jupyter-lab","--allow-root","--ip=0.0.0.0","--port=8888","--no-browser","--NotebookApp.token=''","--notebook-dir=/workspace"]