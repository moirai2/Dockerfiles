FROM moirai2/r
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
samtools \
bedtools \
perl \
pandoc \
&& R -e "install.packages('dplyr',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('Rmisc',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggplot2',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('knitr',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('magicfor',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('reshape2',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('tclust',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('rmarkdown',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
