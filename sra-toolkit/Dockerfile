FROM ncbi/sra-tools
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ENV DEBIAN_FRONTEND noninteractive
ADD vdb-config.txt /root
RUN vdb-config -i --interactive-mode textual << /root/vdb-config.txt \
&& rm /root/vdb-config.txt
