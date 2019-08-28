library(ggplot2)


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

ggplot(data, aes(x = V3, y = V4, colour = V5)) +
  geom_point() +
  facet_wrap( ~ V5)

ggplot(sub_data2, aes(x = V3, y = V4, colour = V5)) +
  geom_point() 
+
  facet_wrap( ~ V5)
