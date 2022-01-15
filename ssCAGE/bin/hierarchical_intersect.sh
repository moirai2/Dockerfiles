#!/bin/sh

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

if [ $# -lt 2 ]
then
echo "hierarchical_intersect.sh '[option]' [input_bed] [out_dir] [annotation_bed] [annotation_bed2] .."
exit
fi

option=$1; shift;
input_file=$1; shift;
output_directory=$1; shift;
if [ "${output_directory}" = "" ]
then
output_directory="output"
fi
if [ ! -e ${output_directory} ]
then
mkdir ${output_directory}
fi

input_bed=${input_file}
input_basename=`basename $input_bed .bed`
processed_once=0
for current_bed in $@
do
remaining_bed=`mktemp mergeXXXXXXXX`
basename=`basename $current_bed .bed`
basename=`basename ${basename} .gff`
basename=`basename ${basename} .vcf`
echo ">intersectBed -a ${input_bed} -b ${current_bed} -u ${option} > ${output_directory}/${input_basename}.${basename}"
intersectBed -a ${input_bed} -b ${current_bed} -u ${option} > ${output_directory}/${input_basename}.${basename}.bed
echo ">intersectBed -a ${input_bed} -b ${current_bed} -v ${option} > ${remaining_bed}"
intersectBed -a ${input_bed} -b ${current_bed} -v ${option} > ${remaining_bed}
if [ "${processed_once}" = "1" ]
then
echo ">rm ${input_bed}"
rm ${input_bed}
fi
input_bed=${remaining_bed}
processed_once=1
done
echo ">mv ${remaining_bed} ${output_directory}/${input_basename}.intergenic.bed"
mv ${remaining_bed} ${output_directory}/${input_basename}.intergenic.bed
