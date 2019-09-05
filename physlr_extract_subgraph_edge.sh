






























dir="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/"
file_name="f1chr4.edge_weights.labeled_f.tsv"
cat f1chr4.edge_weights.labeled_f.tsv | awk '{if($4>80 && $7>30000 && $9>200000) print}' > ${dir}f1chr4_subset_false_edges.tsv
head ${dir}f1chr4_subset_false_edges.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr subgraph -vCGGAACCTCGCCTGAG-1,CTCTTTCCATCACGAT-1 -d1 f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv 2> subgraph.out&
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=50 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv 2> gv.out &
