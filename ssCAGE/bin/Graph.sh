#!/bin/bash
# Graph.sh [arguments]
# A wrapper script to run my Graph program.
# This is needed because Graph.jar can't be searched through path environmental setting.
# Added get_value_from_moirai_config.pl to get config from Moirai.config
# Akira Hasegawa
# 2012/06/20

# Copyright (C) 2013 Akira Hasegawa <ah3q@gsc.riken.jp>
#
# This file is part of MOIRAI.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# "There is insufficient memory for the Java Runtime Environment to continue.  pthread_getattr_np"
# jps -v
# 6716 Jps -Dapplication.home=/usr/lib/jvm/java-1.6.0-openjdk-1.6.0.0.x86_64 -Xms8m
# There is a need to increase memory allocation with -Xms8g and -Xmx8g

if [ $# -lt 1 ]
then
echo "Graph.sh type options < input > output"
exit
fi

directory=`dirname $0`
jre_path=`$directory/get_value_from_moirai_config.pl 'jre path' < $directory/Moirai.config`
$jre_path -Xms8g -Xmx8g -Djava.awt.headless=true -jar $directory/Graph.jar $@
