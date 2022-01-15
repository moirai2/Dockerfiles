FROM r-base
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
bedtools \
ruby \
parallel \
r-cran-fastica \
r-cran-tidyr \
gawk \
&& apt-get remove -y mawk \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
ADD bin /usr/local/bin
