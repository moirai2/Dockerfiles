FROM moirai2/python3
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN pip install jupyter \
&& rm -rf /root/.cache/ \
&& mkdir -p /root/.jupyter/ \
&& mkdir -p /home/jupyter/
ADD jupyter_notebook_config.py /root/.jupyter/jupyter_notebook_config.py
CMD ["jupyter","notebook"]
