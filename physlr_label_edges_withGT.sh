#!/usr/bin/bash

# Description:
# label edges of an input overlap graph based on the ground truth
# 
if [ -n "$1" ];then
	file_overlap=$1	
else
	file_overlap="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.alledges"
fi
if [ -n "$2" ];then
	file_true_edges=$2		
else
	file_true_edges="/projects/btl_scratch/aafshinfard/projects/physlr2/extra/ground_truth/f1chr4_corrected2_GT.txt"
fi
printf "\n ##################################################################################"
printf "\n #### label edges of an input overlap graph based on the ground truth "
printf "\n Overlap file:\n ${file_overlap}.tsv"
printf "\n Ground truth file:\n ${file_true_edges}"

echo ${file_overlap}

#tail ${file_overlap}.tsv
#tail ${file_true_edges}

awk '(NF>3) {print} {s=NF}' ${file_overlap}.tsv > ${file_overlap}.edges.tsv.

#head ${file_overlap}.edges.tsv
# Labeling
awk 'BEGIN {ORS=""} NR==1{{print} {print "\tlabel\n"}}' ${file_overlap}.edges.tsv > ${file_overlap}.edges.labeled_t.tsv
awk 'BEGIN {ORS=""} NR==1{{print} {print "\tlabel\n"}}' ${file_overlap}.edges.tsv > ${file_overlap}.edges.labeled_f.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$1$2]++;next} ($1$2 in a){{print} {print "\t1\n"}}' ${file_true_edges} ${file_overlap}.edges.tsv >> ${file_overlap}.edges.labeled_t.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$2$1]++;next} ($1$2 in a){{print} {print "\t1\n"}}' ${file_true_edges} ${file_overlap}.edges.tsv >> ${file_overlap}.edges.labeled_t.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$1$2]++;next} !($1$2 in a){{print} {print "\t0\n"}}' ${file_overlap}.edges.labeled_t.tsv ${file_overlap}.edges.tsv >> ${file_overlap}.edges.labeled_f.tsv
#awk 'BEGIN {ORS=""} NR==FNR{a[$2$1]++;next} !($1$2 in a){{print} {print "\t0\n"}}' ${file_overlap}.edges.labeled_t.tsv ${file_overlap}.edges.tsv >> ${file_overlap}.edges.labeled_f.tsv




rm ${file_overlap}.edges.tsv
