

cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}' cat hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | awk '{s+=NF==2}END{print s}'

head -n3462510 hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv

grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv | wc -l
grep '_10\s' vertices_hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.tsv > molecules_10_raw.tsv
cat molecules_10_raw.tsv | awk '{print $1}' > molecules_10_raw2.tsv
mv  molecules_10_raw2.tsv  molecules_10_raw.tsv
cat molecules_10_raw.tsv | sed -e "s/$*_10$//" > molecules_10_barcode.tsv


make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs d=10 v="$(paste -d, -s molecules_10_raw.tsv)"

make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.subgraphs/....gv min_component_size=0 -n (use it in the follwing commands)
cd hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.subgraphs
for i in *.tsv; do env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr filter --min-component-size=0 -Ogv -V1 $i > "$i.gv"; done
for i in *.tsv.gv; do /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Teps -o "$i.sfdp.eps" $i; done
for i in *sfdp.eps; do convert $i "$i.pdf"; done

subdir=hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/
command time -v -o "${subdir}1212_446_74.backbobe.map.grch38.n10.paf.gz.time"  env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr map-paf -V$V -n10 "${subdir}1212_446_74.backbobe.path" hg004.k40-w32.n100-5000.c2-x.physlr.tsv grch38/grch38.k40-w32.physlr.tsv | pigz -p16 > "${subdir}1212_446_74.backbobe.map.grch38.n10.paf.gz"

# make backbone.tsv
for i in ${subdir}*[0-9].tsv; do command time -v -o "${i}.backbone.tsv.time" env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr backbone-graph --prune-branches=0 --prune-bridges=0 --prune-junctions=0 -s0 -V1 "${i}" >"${i}.backbone.tsv"; done

for i in *.backbone.tsv; do backbone.path
for i in *.backbone.tsv; do command time -v -o "${subdir}${i}.map.grch38.n10.paf.gz.time"  env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr map-paf -V1 -n10 "${subdir}${i}.path" hg004.k40-w32.n100-5000.c2-x.physlr.tsv grch38/grch38.k40-w32.physlr.tsv | pigz -p16 > "${subdir}${i}.map.grch38.n10.paf.gz"
