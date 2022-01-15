FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
ca-certificates \
curl \
&& cd /root \
&& curl -sSL -o output.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
&& tar -zxf output.tar.gz \
&& rm output.tar.gz \
&& mv sratoolkit.*-ubuntu64/bin/* /usr/local/bin/. \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*