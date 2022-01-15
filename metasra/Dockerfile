FROM moirai2/python3
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN ln -fs /usr/share/zoneinfo/Asia/Tokyo /etc/localtime \
&& apt-get update \
&& apt-get install -y --no-install-recommends \
g++ \
&& pip install \
numpy \
scipy \
scikit-learn \
setuptools \
marisa-trie \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& rm -rf /root/.cache/
