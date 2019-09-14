






























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
# /projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# so:
# cp /projects/btl/jowong/github/physlr/data_n_with_w_filtering/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv
# then:

tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv		#n10
tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv 	#n10 w15000

## n10 w0: all edges
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=50 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.gv 2> gv.out&
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.gv > gv.sfdp.out&
## n10 w15000: all edges
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=50 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.gv 2> gv.n0.out&
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.gv > gv.n0.sfdp.out&

## n10 w0: true edges only
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=50 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv 2> gv.true_edges.out&
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.true_edges.gv > gv.sfdp.true_edges.out&

## n10 w15000: true edges only

awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
awk 'NR==FNR{c[$2,$1]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=50 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.gv 2> gv.n0.true_edges.out&
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.gv > gv.n0.sfdp.true_edges.out&

## n10 and no w15000
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
awk '(s>3 && NF>3 && $4>15000) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv
make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.gv -n

### subset edges - w15000
# 50 %
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.50pe.tsv
awk '(s>2 && NF>2 && NR%2>0) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.50pe.tsv
file_graph="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.50pe"
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv ${file_graph}.tsv >${file_graph}.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o ${file_graph}.gv.sfdp.pdf ${file_graph}.gv > sfdp.out&
# 25 %
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.25pe.tsv
awk '(s>2 && NF>2 && NR%4>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.25pe.tsv
file_graph="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.25pe"
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv ${file_graph}.tsv >${file_graph}.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o ${file_graph}.gv.sfdp.pdf ${file_graph}.gv > sfdp.out&
# 10 %
file_graph="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.10pe"
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv > ${file_graph}.tsv
awk '(s>2 && NF>2 && NR%10>8) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.w15000.tsv >> ${file_graph}.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv ${file_graph}.tsv >${file_graph}.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o ${file_graph}.gv.sfdp.pdf ${file_graph}.gv > sfdp.out&



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



### Making a labeled dataset
#
file_overlap="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap"
file_true_edges="/projects/btl/jowong/github/physlr/ground_truth/true_edges.txt"
tail ${file_overlap}.tsv
tail ${file_true_edges}

awk '(NF>3 && s>3) {print} {s=NF}' ${file_overlap}.tsv > ${file_overlap}.edges.tsv
head ${file_overlap}.edges.tsv
# Labeling
awk 'NR==FNR{a[$1$2]++;next} ($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t1"}' ${file_true_edges} ${file_overlap}.edges.tsv > ${file_overlap}.edges.labeled_t.tsv
awk 'NR==FNR{a[$2$1]++;next} ($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t1"}' ${file_true_edges} ${file_overlap}.edges.tsv >> ${file_overlap}.edges.labeled_t.tsv
awk 'NR==FNR{a[$1$2]++;next} !($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t0"}' ${file_overlap}.edges.labeled_t.tsv ${file_overlap}.edges.tsv > ${file_overlap}.edges.labeled_f.tsv

### With 5 columns
#awk 'NR==FNR{a[$1$2]++;next} ($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t1"}' ${file_true_edges} f1chr4.nn.overlap.edges.tsv > f1chr4.nn.overlap.edges.labeled_t.tsv
#awk 'NR==FNR{a[$2$1]++;next} ($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t1"}' ${file_true_edges} f1chr4.nn.overlap.edges.tsv >> f1chr4.nn.overlap.edges.labeled_t.tsv
#awk 'NR==FNR{a[$1$2]++;next} !($1$2 in a){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t0"}' f1chr4.nn.overlap.edges.labeled_t.tsv f1chr4.nn.overlap.edges.tsv > f1chr4.nn.overlap.edges.labeled_f.tsv


### Making a labeled dataset 3
#
file_weights="/projects/btl/jowong/github/physlr/data_n_with_w_filtering_test/histogramall.txt"
file_true_edges="/projects/btl/jowong/github/physlr/ground_truth/true_edges.txt"
file_output_prefix="f1chr4.edge_weights"
tail ${file_weights}
tail ${file_true_edges}

# Labeling
cat ${file_weights} -n1 | awk 'BEGIN {ORS=""} {print; print"\tlabel\n"}' > ${file_output_prefix}.labeled_t.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$1$2]++;next} ($1$2 in a){print; print "\t1\n"}' ${file_true_edges} ${file_weights} >> ${file_output_prefix}.labeled_t.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$2$1]++;next} ($1$2 in a){print; print "\t1\n"}' ${file_true_edges} ${file_weights} >> ${file_output_prefix}.labeled_t.tsv
cat ${file_weights} -n1 | awk 'BEGIN {ORS=""} {print; print"\tlabel\n"}' > ${file_output_prefix}.labeled_f.tsv
awk 'BEGIN {ORS=""} NR==FNR{a[$1$2]++;next} !($1$2 in a){print; print "\t0\n"}' ${file_output_prefix}.labeled_t.tsv ${file_weights} >> ${file_output_prefix}.labeled_f.tsv
#filtering for n<10
cat ${file_output_prefix}.labeled_f.tsv | awk '(NR<2 || ( NR>1 && $3>10 ) ){print}' > ${file_output_prefix}.labeled_f2.tsv
mv ${file_output_prefix}.labeled_f2.tsv ${file_output_prefix}.labeled_f.tsv





### removing false edges from subgraph ".tsv" files
file_subgraph="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph_CGGAACCTCGCCTGAG-1,CATTCTACAGGGCATA-1"
file_true_edges="/projects/btl/jowong/github/physlr/ground_truth/true_edges.txt"
awk '(NF<3) {print}' ${file_subgraph}.tsv > ${file_subgraph}_t.tsv
awk '(NF>3 && s<3) {print} {s=NF}' ${file_subgraph}.tsv >> ${file_subgraph}_t.tsv
awk '(NF>3 && s>3) {print} {s=NF}' ${file_subgraph}.tsv > ${file_subgraph}.edges.tsv
awk 'NR==FNR{a[$1$2]++;next} ($1$2 in a){print}' ${file_true_edges} ${file_subgraph}.edges.tsv >> ${file_subgraph}_t.tsv
awk 'NR==FNR{a[$2$1]++;next} ($1$2 in a){print}' ${file_true_edges} ${file_subgraph}.edges.tsv >> ${file_subgraph}_t.tsv


