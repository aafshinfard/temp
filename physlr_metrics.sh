for n in $(seq 167);
do
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.paf.gz.pdf  k=40 w=32 n=${n} min_component_size=200   ref=grch38 lr=hg004 r=/projects/btl/jowong/github/physlr/data_hg004_ntHits/ntHits/hg004_k40_ntHits_depth60_repetitive_kmer.bf
make hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.chain.paf.gz  k=40 w=32 n=${n} min_component_size=200   ref=grch38 lr=hg004 r=/projects/btl/jowong/github/physlr/data_hg004_ntHits/ntHits/hg004_k40_ntHits_depth60_repetitive_kmer.bf

pigz -p16 -cd hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.chain.paf.gz | env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl/jowong/github/physlr/bin/physlr liftover_paf grch38/grch38.k40-w32.physlr_pos.tsv - | pigz -p 16 > hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.post_lift_over.chain.paf.gz

pigz -p16 -cd hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.post_lift_over.chain.paf.gz | env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr paf_metrics -g 3000000000 - > hg004.k40-w32.n100-5000.c2-x.physlr.overlap.n${n}.k3.mol.backbone.map.grch38.n10.chain.paf.metrics
done;




10
3
5
3
2
3
3
2
2
3
5
2
2
3
3
2
2
4
3
2
6

