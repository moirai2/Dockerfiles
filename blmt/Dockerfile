FROM moirai2/c
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
ADD blmt_version2 /root/src
RUN cd /root/src \
&& make \
&& mv /root/src/faa2srt /usr/local/bin/. \
&& mv /root/src/ngrams /usr/local/bin/. \
&& mv /root/src/proteinCount /usr/local/bin/. \
&& mv /root/src/proteinNGram /usr/local/bin/. \
&& mv /root/src/srt2lcp /usr/local/bin/. \
&& rm -r /root/src