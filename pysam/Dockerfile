FROM moirai2/python3
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
&& apt-get update \
&& apt-get install -y --no-install-recommends \
libbz2-dev \
zlib1g-dev \
liblzma-dev \
libreadline-dev \
libssl-dev \
libncurses5-dev \
xz-utils \
tk-dev \
libffi-dev \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& pip install cython \
&& pip install wheel \
&& pip install pysam \
&& rm -rf /root/.cache/
