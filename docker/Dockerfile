FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
curl \
build-essential \
&& rm -rf /var/lib/apt/lists/* \
&& curl -ksSL https://get.docker.com/ | sh
