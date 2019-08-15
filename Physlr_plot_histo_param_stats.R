


# first extract the edge weights from the overlap graph file you have (*.overlap.tsv) via the command pattern below:
# cat hg004.k32-w32.n100-5000.c2-x.physlr.overlap.tsv | awk '{if(gsub(/\t/,"\t")>1) print $3}' > hg004.k32-w32.n100-5000.c2-x.physlr.overlap.edge_weights
# then 
library("ggplot2")
data<-read.csv("/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.edge_weights", sep = "\t" , header=T)
head(data)
class(data[1,1])
hist(data[,1])
dim(data)
df = data.frame(data)
head(df)
sub_df = data.frame(df[sample(nrow(df), 1000000), ])
class(sub_df)
names(sub_df) = "edge_weight"
head(sub_df)
p = ggplot(sub_df, aes(x=edge_weight))+geom_histogram(binwidth = 1)+xlim(0,100)
