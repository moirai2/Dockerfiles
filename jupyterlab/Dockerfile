FROM python:latest
USER root
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV TZ JST-9
RUN apt-get update \
&& apt-get -y install \
vim \
less \
&& apt-get autoremove -y \
&& apt-get clean \
&& rm -rf /usr/local/src/* \
&& mkdir -p /root/src
COPY requirements.txt /root/src
WORKDIR /root/src
RUN pip install --default-timeout=100 --upgrade --no-cache-dir \
pip \
black \
import-ipynb \
ipywidgets \
jupyterlab \
jupyterlab-git \
jupyterlab_code_formatter \
jupyterlab_widgets \
lckr-jupyterlab-variableinspector \
setuptools \
&& rm -rf ~/.cache/pip
RUN pip install --default-timeout=100 --upgrade --no-cache-dir -r requirements.txt \
&& rm -rf ~/.cache/pip
RUN pip install --default-timeout=100 --upgrade --no-cache-dir bash_kernel \
&& python -m bash_kernel.install
CMD [" jupyter-lab","--allow-root","--ip=0.0.0.0","--port=8888","--no-browser","--NotebookApp.token=''","--notebook-dir=/workspace"]