library(ggplot2)
#############################################################
### plot 2 features (n and w)

ad_t = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.edges.labeled_t.tsv"
ad_f = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.edges.labeled_f.tsv"
data_t = read.table(ad_t, sep = "\t", header = FALSE)
data_f = read.table(ad_f, sep = "\t", header = FALSE)
# subsample 
sub_t = data.frame(data_t[sample(nrow(data_t), dim(data_f)[1]), ])
data = rbind(sub_t,data_f)

dim(data_t)
dim(data)
head(data)
tail(data)
class(data[1,3])
data$V5 <- as.factor(data$V5)

sub_data2 = data.frame(data[sample(nrow(data), 30000), ])

ggplot(data, aes(x = n, y = w, colour = "true edges")) +
  geom_point() +
  facet_wrap( ~ V5)
names(sub_data2) = c("V1", "V2", "V3", "V4", "V5")
ggplot(sub_data2, aes(x = V3, y = V4, colour = V5)) +
  geom_point() + ggtitle("f1chr4 n vs. w") + xlab("n") + ylab("w")+ theme(text = element_text(size=16))
+ facet_wrap( ~ V5)

#############################################################
### plot all features included in 2 files.

ad_t = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.edge_weights.labeled_t.tsv"
ad_f = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.edge_weights.labeled_t.tsv"
data_t = read.table(ad_t, sep = "\t", header = TRUE)
data_f = read.table(ad_f, sep = "\t", header = TRUE)
# subsample 
sub_t = data.frame(data_t[sample(nrow(data_t), dim(data_f)[1]), ])
data = rbind(sub_t,data_f)

dim(data_t)
dim(data)
head(data)
tail(data)
class(data[1,3])
data$V5 <- as.factor(data$V5)

sub_data2 = data.frame(data[sample(nrow(data), 30000), ])

ggplot(data, aes(x = V3, y = V4, colour = V5)) +
  geom_point() +
  facet_wrap( ~ V5)

ggplot(sub_data2, aes(x = V3, y = V4, colour = V5)) +
  geom_point() 
+
  facet_wrap( ~ V5)

