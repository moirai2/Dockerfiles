FROM ubuntu
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
perl \
libio-socket-ssl-perl \
libhtml-html5-entities-perl \
libhtml-entities-numbered-perl \
libwww-perl \
libxml-simple-perl \
curl \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& echo y | sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
RUN cd /root/edirect/ \
&& curl --output xtract.Linux.gz ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz \
&& gunzip -f xtract.Linux.gz \
&& chmod +x xtract.Linux \
&& tar -zxf Mozilla-CA.tar.gz \
&& mv * /usr/local/bin/.