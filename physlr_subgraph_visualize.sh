

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
for i in *.backbone.tsv; do echo "$i" | j=$(sed -e "s/$.backbone.tsv$//"); echo j
for i in *.backbone.tsv; do command time -v -o "${subdir}${i}.map.grch38.n10.paf.gz.time"  env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr map-paf -V1 -n10 "${subdir}${i}.path" hg004.k40-w32.n100-5000.c2-x.physlr.tsv grch38/grch38.k40-w32.physlr.tsv | pigz -p16 > "${subdir}${i}.map.grch38.n10.paf.gz"

gunzip -c hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/1212_446_74.backbone.map.grch38.n10.paf.gz | env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr pypy3 /projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr annotate-graph -V1 --min-component-size=0 -Ogv -V1 hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/1212_446_74.backbone.tsv hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/1212_446_74.backbone.tsv - >hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/1212_446_74.backbone.map.grch38.n10.ann.gv



###########################

pattern: for i in *tsv; do echo $i | sed -e "s/$*.tsv$//" | head -n1 "$(cat).tsv"; done;
### variables:
# Path to the Physlr project.
physlr_path="/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/"
#physlr_path=$(pwd)
# Path to the directory for profiling
path_profile=$(PWD)
# Path to the Physlr executable.
bin=${physlr_path}bin
python_executable=/projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3
# Python interpreter.
python="env PYTHONPATH=${physlr_path} ${python_executable}"
subir=hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs/


# Find barcodes/molecules of interest

# Make subgraphs using barcodes/molecules
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.subgraphs d=10 v="$(paste -d, -s molecules_10_barcode.tsv)"
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.subgraphs d=10 v="$(paste -d, -s molecules_10_raw.tsv)"

# fix naming (.mol)
for i in ${subdir}*[0-9].tsv; do echo $i | sed -e "s/$*.tsv$//" | cp ${i} "$(cat).mol.tsv"; done

# make gv file from tsv file
for i in ${subdir}*[0-9].tsv; do echo $i | sed -e "s/$*.tsv$//" | $python ${bin}/physlr filter --min-component-size=0 -Ogv -V1 $i > "${subdir}$(cat).gv"; done
#/projects/btl_scratch/aafshinfard/projects/physlr/publication/jowong/hg004/physlr/bin/physlr

#make backbone file from tsv file
for i in ${subdir}*[0-9].mol.tsv; do echo $i | sed -e "s/$*.tsv$//" | $python ${bin}/physlr backbone-graph --prune-branches=10 --prune-bridges=0 --prune-junctions=0 -s0 -V1 ${i} >"$(cat).backbone.tsv"; done
for i in ${subdir}*[0-9].mol.backbone.tsv; do echo $i | sed -e "s/$*.tsv$//" | $python ${bin}/physlr backbone --prune-branches=0 -s0 -V1 ${i} >"$(cat).path"; done

# Map
for i in ${subdir}*[0-9].backbone.path; do echo $i | sed -e "s/$*.backbone.path$//" | $python ${bin}/physlr map-paf -V1 -n10 "$(cat).backbone.path" hg004.k40-w32.n100-5000.c2-x.physlr.overlap.m85.mol.split.tsv grch38/grch38.k40-w32.physlr.tsv | pigz -p16 > "${i}.map-split.grch38.n10.paf.gz"; done

# Map the reference to the backbone graph and output PAF.
%.map.$(ref).n10.paf.gz: %.path $(lr).k$k-w$w.n$(minimum_barcode_multiplicity)-$(maximum_barcode_multiplicity).c2-x.physlr.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf -V$V -n10 $^ | $(gzip) >$@

# Map the reference to the backbone graph with split minimizers and output PAF.
%.backbone.map-split.$(ref).n10.paf.gz: %.backbone.path %.split.tsv $(name)/$(ref).k$k-w$w.physlr.tsv
	$(time) $(python) $(bin)/physlr map-paf --mx-type split -V$V -n10 $^ | $(gzip) >$@
  
  
for i in ${subdir}*.backbone.tsv.map-split.grch38.n10.paf.gz; do echo $i | sed -e "s/$*.backbone.tsv.map-split.grch38.n10.paf.gz$//" | $python ${bin}/physlr annotate-graph -V$V --min-component-size=$(min_component_size) -Ogv -V$V $*.backbone.tsv $*.backbone.path - >$@

%.backbone.map.$(ref).n10.ann.gv: %.backbone.map.$(ref).n10.paf.gz %.backbone.tsv %.backbone.path
	gunzip -c $< | $(python) $(bin)/physlr annotate-graph -V$V --min-component-size=$(min_component_size) -Ogv -V$V $*.backbone.tsv $*.backbone.path - >$@
