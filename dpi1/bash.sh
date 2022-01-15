#!/bin/sh

user=`whoami`
workdir=`pwd`
docker run -it --rm -v $workdir:/home/$user --workdir /home/$user moirai2/dpi1 bash
