FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
python3 \
python3-pip \
python3-setuptools \
python3-dev \
gcc \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& ln -s /usr/bin/python3 /usr/bin/python \
&& ln -s /usr/bin/pip3 /usr/bin/pip
