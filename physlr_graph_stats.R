
library("ggplot2")
library("igraph")


dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n25.3.stats"
dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n40.2.stats"
dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n70.3.stats"

# s_t = -1 + (25 - 10) / 2
a = read.table(dat, header = TRUE, fill = TRUE)
# a[1:s_t,"Large.comp.count"]=0
# a[1:s_t,"Comp.count"]=0
# a[1:s_t,"Size.largest"]=0


head(a)
tail(a)
ggplot(data = a, aes(x=Threshold, y = Comp.count)) + geom_col(fill='purple', color = 'purple') + labs(x="Threshold on edge weight", y="Number of components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Large.comp.count)) + geom_col(fill='purple', color = 'purple') + labs(x="Threshold on edge weight", y="Number of significant components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Size.largest)) + geom_col(fill='purple', color = 'purple') + labs(x="Threshold on edge weight", y="Size of largest component (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Edge.count)) + geom_col(fill='purple', color = 'purple') + labs(x="Threshold on edge weight", y="Number of edges (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Node.count)) + geom_col(fill='purple', color = 'purple') + labs(x="Threshold on edge weight", y="Number of vertices (in the overlap graph)")


ggplot(data = a, aes(x=Threshold, y = Node.count)) + geom_col()
ggplot(data = a, aes(x=Threshold, y = Edge.count)) + geom_col()

