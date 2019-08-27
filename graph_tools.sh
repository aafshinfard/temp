






























file="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.2.tsv"
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.tsv | awk '{ if (NF>3) {print $1"\t"$2"\t"$3} else if (NF==2) {print $1"\t"$2} else if (1) print}' > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.tsv
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.tsv | awk '{ if (NF>3) {print $1"\t"$2"\t"$4} else if (NF==2) {print $1"\t"$2} else if (1) print}' > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.tsv
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.tsv | sed s/w/n/ > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.2.tsv

nohup make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.gv &
nohup make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.2.gv &

cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.tsv | awk 'NF==0 {print s; print; getline; print} {s=$1$2}'

#ccomps connected components of the graph
ccomps -xz f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.gv -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.ccomps.gv


###
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n51.mol.tsv | sed s/_0/""/g > cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n51.edited.mol.tsv


# awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n51.mol.edited.tsv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.tsv > mixed.txt
f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.tsv


# file with no W filtering:
/projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# so:
cp /projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# then:
## n10 and no W filtering + true edges only
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv -n 
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv.sfdp.pdf -n

## n10 and no w15000
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
awk '(s>3 && NF>3 && $4>15000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.gv -n

## n10 and no w20000
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
awk '(s>3 && NF>3 && $4>17000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.gv -n
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.gv.sfdp.pdf -n


awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv -n 
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv.sfdp.pdf -n



tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv
## 

##############################################################################
### Fly overlap graphs

# file with no W filtering:
/projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# so:
cp /projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# then:







tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv		#n10
tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv 	#n10 w15000

## n10 w0: all edges
f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
## n10 w15000: all edges
f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv

## n10 w0: true edges only
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv -n 
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv.sfdp.pdf -n

## n10 w15000: true edges only

awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.gv -n 
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.gv.sfdp.pdf -n






## n10 and no w15000
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
awk '(s>3 && NF>3 && $4>15000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.gv -n

## n10 and no w17000
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
awk '(s>3 && NF>3 && $4>17000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.gv -n
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w17000.gv.sfdp.pdf -n


awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv -n 
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv.sfdp.pdf -n



tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv
## 

awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv


awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.w15000.tsv
awk '(s>3 && NF>3 && $4>15000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.w15000.tsv

awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv


awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.w15000.tsv
awk '(s>3 && NF>3 && $4>15000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.w15000.tsv

