#!/usr/bin/bash

file_w=$1
file_wj=$2
threshold=$3

cat $1 | awk '{if($3 >= 15000) print $1,$2,$3}' > w_above_${threshold}.txt
cat $1 | awk '{if($3 < 15000) print $1,$2,$3}' > w_below_${threshold}.txt

awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' w_above_${threshold}.txt $2 > wj_above_${threshold}.txt
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' w_below_${threshold}.txt $2 > wj_below_${threshold}.txt
