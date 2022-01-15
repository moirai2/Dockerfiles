FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
perl \
build-essential \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
