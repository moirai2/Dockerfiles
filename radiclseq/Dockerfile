FROM moirai2/r
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y \
samtools \
bedtools \
libcurl4-gnutls-dev \
libxml2-dev \
libssl-dev \
&& R -e "install.packages('AnnotationHub',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('BSgenome.Mmusculus.UCSC.mm10',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('Cairo',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('DT',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('GGally',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('GenomicFeatures',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('RColorBrewer',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('biomaRt',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('cowplot',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('data.table',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggforce',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggplot2',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggpubr',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggrepel',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggsci',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('ggthemes',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('gridExtra',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('knitr',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('png',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('scales',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('viridis',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& R -e "install.packages('rtracklayer',dependencies=TRUE,repos='http://cran.rstudio.com/')" \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/*
