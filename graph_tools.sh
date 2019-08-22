






























file="f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.2.tsv"
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.tsv | awk '{ if (NF>3) {print $1"\t"$2"\t"$3} else if (NF==2) {print $1"\t"$2} else if (1) print}' > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.tsv
cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.tsv | awk '{ if (NF>3) {print $1"\t"$2"\t"$4} else if (NF==2) {print $1"\t"$2} else if (1) print}' > f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.tsv


nohup make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisW.gv &
nohup make f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.gv &

cat f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.n15000.nisn.tsv | awk 'NF==0 {print s; print; getline; print} {s=$1$2}'

ccomps #for connected components of the graph
ccomps -xz 
