FROM ubuntu
LABEL  maintainer "Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
perl \
curl \
rna-star \
samtools \
bedtools \
fastqc \
kraken \
libio-socket-ssl-perl \
libhtml-html5-entities-perl \
libhtml-entities-numbered-perl \
libwww-perl \
libxml-simple-perl \
sra-toolkit \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& echo y | sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
