

cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}' cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}'

head -n3462510 hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv

grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | wc -l
grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > molecules_10_raw.tsv
cat molecules_10_raw.tsv | awk '{print $1}' > molecules_10_raw2.tsv
mv  molecules_10_raw2.tsv  molecules_10_raw.tsv
cat molecules_10_raw.tsv | sed -e "s/$*_10$//" > molecules_10_barcode.tsv

***exclude_source=1***
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs d=10 v="$(paste -d, -s molecules_10_raw.tsv)"




