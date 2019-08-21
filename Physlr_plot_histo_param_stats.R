


# first extract the edge weights from the overlap graph file you have (*.overlap.tsv) via the command pattern below:
# cat hg004.k32-w32.n100-5000.c2-x.physlr.overlap.tsv | awk '{if(gsub(/\t/,"\t")>1) print $3}' > hg004.k32-w32.n100-5000.c2-x.physlr.overlap.edge_weights
# then 
library("ggplot2")
data<-read.csv("/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.k32-w32.n100-5000.c2-x.physlr.overlap.edge_weights", sep = "\t" , header=T)
#Jonathans output:
data<-read.csv("/projects/btl/jowong/github/physlr/data/histogram.txt", sep = "\t" , header=T)
data<-read.csv("/projects/btl/jowong/github/physlr/data/histogramv4.txt", sep = "\t" , header=T)
data<-read.csv("/projects/btl/jowong/github/physlr/data_jaccard/histogramj.txt", sep = "\t" , header=T)


data<-read.csv("/projects/btl/jowong/github/physlr/data_n_default10/histogramw.txt", sep = "\t" , header=T)
data_backup = data
head(data)
# if you data hsa 3 columns, discard the first two using:
data= data.frame(data[,3])





class(data[1,1])
#hist(data[,1])
dim(data)
df = data.frame(data)
head(df)
sub_df = data.frame(df[sample(nrow(df), 1000000), ])
class(sub_df)
names(sub_df) = "edge_weight"
head(sub_df)
p = ggplot(sub_df, aes(x=edge_weight))+geom_histogram(binwidth = 100)#+xlim(0,100)
p2 = ggplot(sub_df, aes(x=edge_weight))+geom_histogram(binwidth = 100)+xlim(0,50000)+ylim(0,7000)
#p10 = ggplot(round(sub_df/10), aes(x=edge_weight))+geom_histogram(binwidth = 1)#+xlim(0,100)



# Mixed:
data2<-read.csv("/projects/btl/aafshinfard/projects/tempOnGithub/temp/wj_above_15000.txt", sep = "\t" , header=F)
data3<-read.csv("/projects/btl/aafshinfard/projects/tempOnGithub/temp/wj_below_15000.txt", sep = "\t" , header=F)
data_backup2 = data2 # data2 = data_backup2
data_backup3 = data3 # data3 = data_backup3
dim(data2)
dim(data3)
# data2=data3
# if you data hsa 3 columns, discard the first two using:
data2= data.frame(data2[,3])
data3= data.frame(data3[,3])

df2 = data.frame(data2)
df3 = data.frame(data3)
head(df2)
sub_df2 = data.frame(df2[sample(nrow(df2), 1000000), ])
sub_df3 = df3
class(sub_df2)
names(sub_df2) = "edge_weight"
head(sub_df2)
class(sub_df3)
names(sub_df3) = "edge_weight"
head(sub_df3)



p21 = ggplot(sub_df2, aes(x=edge_weight))+geom_histogram(binwidth = 0.005)#+xlim(0,100)
p22 = ggplot(sub_df2, aes(x=edge_weight))+geom_histogram(binwidth = 100)+xlim(0,50000)+ylim(0,7000)

p31 = ggplot(sub_df3, aes(x=edge_weight))+geom_histogram(binwidth = 0.005)#+xlim(0,100)
p32 = ggplot(sub_df3, aes(x=edge_weight))+geom_histogram(binwidth = 100)+xlim(0,50000)+ylim(0,7000)
