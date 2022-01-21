FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
wget \
ca-certificates \
&& wget https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz \
&& tar zxvf express-1.5.1-linux_x86_64.tgz \
&& mv express-1.5.1-linux_x86_64/express /usr/local/bin/. \
&& rm express-1.5.1-linux_x86_64.tgz \
&& rm -r express-1.5.1-linux_x86_64 \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
