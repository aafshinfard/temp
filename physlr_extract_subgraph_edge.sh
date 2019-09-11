






























dir="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/"
file_name="f1chr4.edge_weights.labeled_f.tsv"
cat f1chr4.edge_weights.labeled_f.tsv | awk '{if($4>80 && $7>30000 && $9>200000) print}' > ${dir}f1chr4_subset_false_edges.tsv
head ${dir}f1chr4_subset_false_edges.tsv
## make the subgraph of a specific edge
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr subgraph -vCGGAACCTCGCCTGAG-1,CTCTTTCCATCACGAT-1 -d1 f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv 2> subgraph.out&
## overlap.subgraph
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv > sfdp.out&
## overlap.subgraph.n50
/gsc/btl/linuxbrew/bin/mlr --tsvlite filter '$n >= 50' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv > sfdp.out&

U	V	n_orig	n	n_jaccard	n_tfidf	w	w_jaccard	w_tfidf	label
AAAGCGGAGTACCTGT-1	AGGAAGCTCGTCTGCT-1	35	35	0.03791982665222102	182.28773950380614	18409	0.26509511397837077	117606.88632435875	0
AAAGCGGAGTACCTGT-1	TGCACGGGTGAACAAT-1	41	41	0.04904306220095694	216.9925155147764	18259	0.6343454697053919	116835.31704039147	0
AAAGCGGAGTACCTGT-1	GGAATAATCTGCTACC-1	26	26	0.043993231810490696	144.9375838863627	18571	0.44381512283720487	118486.85978864846	0
AAAGCGGAGTACCTGT-1	GCCATCTGTGTTGTGT-1	17	17	0.016472868217054265	89.80071801017229	18564	0.2515958528156129	118425.71889795162	0
AAAGCGGAGTACCTGT-1	TCGCGTTAGTTGTAGA-1	14	14	0.029288702928870293	74.0365767537969	18223	0.6909980282117397	116652.3563694298	0
CGGAACCTCGCCTGAG-1	CTTGTGCAGTTGGCTT-1	113	113	0.07543391188251002	602.8244552544769	23258	0.4371722336046315	151529.0998789517	0
CGGAACCTCGCCTGAG-1	TTAAGGCAGCTAAACA-1	55	55	0.03873239436619718	284.32249633240747	23124	0.4297981487677038	149077.55231691446	0
CGGAACCTCGCCTGAG-1	CTCTTTCCATCACGAT-1	102	102	0.06850235057085292	534.3541574389083	32515	0.4059554279293339	204438.79992588575	0
CGGAACCTCGCCTGAG-1	CCTCATGGTGCTAGCC-1	97	97	0.06764295676429567	502.08556145981294	37352	0.5039735546110774	229608.21134964406	0
CGGAACCTCGCCTGAG-1	CGACGTGCAGTATCTG-1	43	43	0.046236559139784944	223.3286350651438	22184	0.44122677910815866	145271.148554806	0
CGGAACCTCGCCTGAG-1	GACTACAAGCTGCCCA-1	136	136	0.09590973201692525	714.1534901930787	32013	0.6257427677873338	201816.4206417801	0
CGGAACCTCGCCTGAG-1	TCATGTTAGGTACAAT-1	90	90	0.06298110566829951	479.99973202171725	26256	0.4549960142792777	170264.3565175567	0
CGGAACCTCGCCTGAG-1	CATTCTACAGGGCATA-1	142	142	0.09957924263674614	760.603359615834	30805	0.45661389778252104	195266.06545490728	0
## make the subgraph of a specific edge
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr subgraph -vAAAGCGGAGTACCTGT-1,CTCTTTCCATCACGAT-1 -d1 f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv 2> subgraph.out&
## overlap.subgraph
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.gv > sfdp.out&
## overlap.subgraph.n50
/gsc/btl/linuxbrew/bin/mlr --tsvlite filter '$n >= 50' f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.tsv
nohup env PYTHONPATH=/projects/btl_scratch/aafshinfard/projects/physlr2/physlr /projects/btl/aafshinfard/virtuEnv/pypy3/bin/pypy3 /projects/btl_scratch/aafshinfard/projects/physlr2/physlr/bin/physlr filter --min-component-size=1 -Ogv f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.tsv >f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv 2> gv.out &
nohup /gsc/btl/linuxbrew/bin/sfdp -Gsize=100 -Goverlap_scaling=200 -Tpdf -o f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv.sfdp.pdf f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.subgraph.n50.gv > sfdp.out&

