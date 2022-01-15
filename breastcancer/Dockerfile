FROM moirai2/python3
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
&& pip install \
numpy \
pandas \
plotly \
scipy \
matplotlib \
intervals \
&& rm -rf /root/.cache/
