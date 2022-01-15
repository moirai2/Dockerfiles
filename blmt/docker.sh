#!/bin/bash
image=moirai2/blmt
user=`whoami`
workdir=`pwd`;
homedir=/home/$user
docker run \
  -it \
  --rm \
  -v $workdir:$homedir \
  --workdir $homedir \
  $image \
  bash
