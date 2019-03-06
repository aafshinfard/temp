#!/usr/bin/bash

#awk {'if($2>1000) print $1'} f1.n100-2000.physlr.overlap.n20.tsv > b1000.tsv

query=$1
filename=$2
grep $query ${filename}.tsv > ${query}_${filename}.tsv

	



