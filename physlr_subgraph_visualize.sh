

cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}' cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}'

head -n3462510 hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv

grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | wc -l
grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > molecules_10_raw.tsv



