FROM moirai2/r
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
samtools \
bedtools \
rna-star \
perl \
default-jre  \
&& R -e "install.packages('scales',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('dplyr',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggplot2',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('reshape',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
ADD bin /usr/local/bin