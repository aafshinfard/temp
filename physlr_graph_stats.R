
library("ggplot2")
library("igraph")


dat="/projects/btl_scratch/aafshinfard/physlr2/physlr/data/hg004.indexlr.n100-5000.c2-x.physlr.overlap.n10.stats"





a = read.table(dat, header = TRUE, fill = TRUE)
head(a)
ggplot(data = a, aes(x=Threshold, y = Comp.count)) + geom_step()
ggplot(data = a, aes(x=Threshold, y = Node.count)) + geom_col()
ggplot(data = a, aes(x=Threshold, y = Edge.count)) + geom_col()
