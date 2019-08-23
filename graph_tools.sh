






























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



tail f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv
## 
awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.true_edges.tsv


awk 'NF<3 {print} (s<3 && NF>2) {print} {s=NF}' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
awk 'NR==FNR{c[$1,$2]++;next};c[$1,$2]' /projects/btl/jowong/github/physlr/ground_truth/true_edges.txt f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.tsv >> f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
wc -l f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n0.mol.true_edges.tsv
