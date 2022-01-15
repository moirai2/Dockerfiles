FROM php:7.3-apache
LABEL maintainer="Akira Hasegawa <akira.hasegawa@riken.jp>"
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
sqlite3 \
perl \
build-essential \
&& apt -y clean \
&& rm -rf /var/lib/apt/lists/* \
&& curl -sSL https://get.docker.com/ | sh \
&& echo|cpan DBI LWP \
&& rm -rf /root/.cpan \
&& chmod -R 777 /var/www/html
COPY php.ini /usr/local/etc/php/
ADD moirai2.php /var/www/html/
ADD moirai2.pl /var/www/html/
ADD rdf.pl /var/www/html/
ADD js /var/www/html/js/
ADD index.php /var/www/html/
