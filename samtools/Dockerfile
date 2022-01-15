FROM moirai2/c
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libbz2-dev \
zlib1g-dev \
libncurses5-dev \
libncursesw5-dev \
liblzma-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
