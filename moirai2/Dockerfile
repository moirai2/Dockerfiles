# image
FROM php:5.4-apache

# Set php.ini
COPY ./config/php.ini /usr/local/etc/php/
COPY ./config/apache2.conf /etc/apache2/
COPY ./config/sites/*.conf /etc/apache2/sites-available/
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
libpq-dev \
libwww-perl \
perl \
build-essential \
&& docker-php-ext-install pdo_mysql pdo_pgsql mysqli mbstring