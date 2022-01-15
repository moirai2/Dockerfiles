FROM moirai2/jupyter
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN pip3 install wheel \
&& pip3 install tornado \
&& pip install numpy \
&& pip install scikit-learn \
&& pip install matplotlib \
&& pip install seaborn \
&& pip install pandas \
&& pip install umap-learn \
&& rm -rf /root/.cache/
