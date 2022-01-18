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
RUN pip install --default-timeout=100  --upgrade --no-cache-dir \
pip \
setuptools \
jupyterlab \
black \
jupyterlab_code_formatter \
jupyterlab-git \
lckr-jupyterlab-variableinspector \
jupyterlab_widgets \
ipywidgets \
import-ipynb \
&& rm -rf ~/.cache/pip
RUN pip install --default-timeout=100 --upgrade --no-cache-dir -r requirements.txt \
&& rm -rf ~/.cache/pip
CMD [" jupyter-lab","--allow-root","--ip=0.0.0.0","--port=8888","--no-browser","--NotebookApp.token=''","--notebook-dir=/workspace"]



      
      
