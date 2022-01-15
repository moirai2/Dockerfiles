FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
&& apt-get update \
&& apt-get install -y \
apt-transport-https \
software-properties-common \
dirmngr \
&& apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
&& apt-get update \
&& apt-get install -y --no-install-recommends \
build-essential \
r-base \
r-base-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
