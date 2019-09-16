



























from collections import defaultdict
molecule = {}
f = open("/projects/btl/jowong/github/physlr/ground_truth/flchr4.reads.molecule.bed", "r")
#f = open("/projects/btl_scratch/aafshinfard/projects/physlr2/extra/spruce/ws77111.contig2.allreadsbx.molecule.bed", "r")
for i in f:
    columns= i.split("\t")
    molecule[columns[3]] = (int(columns[1]), int(columns[2]))
edges = defaultdict(int)
for i, n in molecule.items():
    for j, m in molecule.items():
        if edges[(j, i)] == 1:
            continue
        if i == j:
            continue
        if (m[0] >= n[0] and m[0]<=n[1]) or (m[1] >= n[0] and m[1]<=n[1]):
            edges[(i,j)] = 1

fout = "edges_GT1.txt"
fo = open(fout, "w")
for k, v in edges.items():
     if v==1:
         fo.write('\n'+str(k[0]) +'\t'+ str(k[1]) +'\t'+str(v))
