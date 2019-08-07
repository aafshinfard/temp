
library("ggplot2")
library("igraph")


dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n25.2.stats"
dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n40.2.stats"
dat="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.n70.3.stats"

a = read.table(dat, header = TRUE, fill = TRUE)
a[1:10,"Large.comp.count"]=0
a[1:10,"Comp.count"]=0
a[1:10,"Size.largest"]=0


head(a)
tail(a)
ggplot(data = a, aes(x=Threshold, y = Comp.count)) + geom_col() + labs(x="Threshold on edge weight", y="Number of components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Large.comp.count)) + geom_col() + labs(x="Threshold on edge weight", y="Number of components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Size.largest)) + geom_col() + labs(x="Threshold on edge weight", y="Number of components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Node.count)) + geom_col()
ggplot(data = a, aes(x=Threshold, y = Edge.count)) + geom_col()

