FROM python
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
&& pip install \
numpy \
pandas \
matplotlib \
&& rm -rf /root/.cache/
