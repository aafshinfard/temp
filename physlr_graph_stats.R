
library("ggplot2")
library("igraph")


dat="/projects/btl_scratch/aafshinfard/physlr2/physlr/data/hg004.indexlr.n100-5000.c2-x.physlr.overlap.n10.2.stats"

a = read.table(dat, header = TRUE, fill = TRUE)
a[1:10,"Comp.count"]=0

head(a)
ggplot(data = a, aes(x=Threshold, y = Comp.count)) + geom_col() + labs(x="n, g, w, or whatever :D", y="Number of components (in the overlap graph)")
ggplot(data = a, aes(x=Threshold, y = Node.count)) + geom_col()
ggplot(data = a, aes(x=Threshold, y = Edge.count)) + geom_col()

