library(ggplot2)
library("Rtsne")
#############################################################
### plot 2 features (n and w)
# first you need to make two files one for true edges and one for false edges and then:
#ad_t = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.edges.labeled_t.tsv"
#ad_f = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.edges.labeled_f.tsv"
ad_t = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.nn.edges.edges.labeled_t.tsv"
ad_f = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.nn.edges.edges.labeled_f.tsv"
data_t = read.table(ad_t, sep = "\t", header = FALSE)
data_f = read.table(ad_f, sep = "\t", header = FALSE)

ad_t ="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.alledges.edges.labeled_t.tsv"
ad_f ="/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.alledges.edges.labeled_f.tsv"
data_t = read.table(ad_t, sep = "\t", header = TRUE)
data_f = read.table(ad_f, sep = "\t", header = TRUE)
head(data_t)
dim(data_t)
head(data_f)
dim(data_f)
hist(data_t[,'n'])
#n10
data_f_n10 = data_f[data_f[,"n"]>=9,]
head(data_f_n10)
dim(data_f_n10)/dim(data_f)
data_t_n10 = data_t[data_t[,"n"]>=9,]
head(data_t_n10)
dim(data_t_n10)/dim(data_t)
(dim(data_f_n10)+dim(data_t_n10))/(dim(data_f)+dim(data_t))
# unbalanced
data_unb = rbind(data_t,data_f)
data_unb_n10 = rbind(data_t_n10,data_f_n10)

# subsample 
sub_t = data.frame(data_t[sample(nrow(data_t), dim(data_f)[1]), ])
sub_t_n10 = data.frame(data_t_n10[sample(nrow(data_t_n10), dim(data_f_n10)[1]), ])
head(sub_t)
data = rbind(sub_t,data_f)
data_n10 = rbind(sub_t_n10,data_f_n10)



dim(data_t)
dim(data)
head(data)
tail(data)
class(data[1,3])
colnames(data)
data$label <- as.factor(data$label)
data_n10$label <- as.factor(data_n10$label)

sub_data2 = data.frame(data[sample(nrow(data), 30000), ])
sub_data2_n10 = data.frame(data_n10[sample(nrow(data_n10), 30000), ])
sub_data2_unb = data.frame(data_unb[sample(nrow(data_unb), 30000), ])
sub_data2_unb_n10 = data.frame(data_unb_n10[sample(nrow(data_unb_n10), 30000), ])
head(sub_data2)

# w histo
ggplot(sub_data2, aes(x = w, colour = label)) +
  stat_bin(bins=30) + ggtitle("f1chr4 w") + xlab("w") + ylab("histo")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(sub_data2_n10, aes(x = w, colour = label)) +
  stat_bin(bins=30) + ggtitle("f1chr4 w") + xlab("w") + ylab("histo")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)


names(sub_data2) = c("V1", "V2", "V3", "V4", "V5")
ggplot(sub_data2, aes(x = V3, y = V5, colour = V6)) +
  geom_point() + ggtitle("f1chr4 n vs. ns1") + xlab("n") + ylab("ns1")+ theme(text = element_text(size=16)) + facet_wrap( ~ V6)
ggplot(sub_data2, aes(x = V3, y = V5)) +
  stat_density_2d()  + ggtitle("f1chr4 n vs. ns1") + xlab("n") + ylab("ns1")+ theme(text = element_text(size=16)) + facet_wrap( ~ V6)

#############################################################
### plot all features included in 2 files.

ad_t = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.edge_weights.labeled_t.tsv"
ad_f = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.edge_weights.labeled_f.tsv"

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
data$label <- as.factor(data$label)

# hitogram of all features

# plot 2 features at a time
sub_data2 = data.frame(data[sample(nrow(data), 30000), ])
head(sub_data2)
ggplot(sub_data2, aes(x = n_tfidf, y = w_tfidf, colour = label)) +
  geom_point() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n_tfidf") + ylab("w_tfidf")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)
ggplot(sub_data2, aes(x = n_tfidf, y = w_tfidf, colour = label)) +
  stat_density_2d() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n_tfidf") + ylab("w_tfidf")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)
ggplot(sub_data2, aes(x = n, y = w, colour = label)) +
  stat_density_2d() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n") + ylab("w")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)

#n10
ggplot(sub_data2_n10, aes(x = n, y = w, colour = label)) +
  stat_density_2d() + ggtitle("f1chr4 n vs. w") + xlab("n") + ylab("w")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)
ggplot(sub_data2_n10, aes(x = n_tfidf, y = w_tfidf, colour = label)) +
  stat_density_2d() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n_tfidf") + ylab("w_tfidf")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)
ggplot(sub_data2_n10, aes(x = n_tfidf, y = w_tfidf, colour = label)) +
  geom_point() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n_tfidf") + ylab("w_tfidf")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)

ggplot(sub_data2_unb_n10, aes(x = n_tfidf, y = w_tfidf, colour = label)) +
  stat_density_2d() + ggtitle("f1chr4 n vs. w (with weighted minimizers)") + xlab("n_tfidf") + ylab("w_tfidf")+ theme(text = element_text(size=16))+ facet_wrap( ~ label)

head(as.matrix(sub_data2[,1:9]))
head((sub_data2[,1:9]))

###############################################
# PCA
pca = prcomp(sub_data2[,3:9], scale=TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var = pca$sdev^2
pca.var.per = round(pca.var/sum(pca.var)*100,1)
barplot(pca.var,main="Scree plot", xlab="PC", ylab="Variation (%)")

pca.data = data.frame(Sample=rownames(pca$x),
                      X=pca$x[,1],
                      Y=pca$x[,2])
head(pca.data)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")
head(pca.data)
head(sub_data2)
pca.data.labeled = cbind(pca.data,sub_data2[,'label'])
names(pca.data.labeled)[4] = "label"
head(pca.data.labeled)
ggplot(data=pca.data.labeled, aes(x=X, y=Y, colour=label))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")


ggplot(pca.data.labeled, aes(x=X, y=Y, colour=label)) +
  geom_point() + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2$label)
ggplot(pca.data.labeled, aes(x=X, y=Y)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2$label)

##### unb:
pca = prcomp(sub_data2_unb[,3:9], scale=TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var = pca$sdev^2
pca.var.per = round(pca.var/sum(pca.var)*100,1)
barplot(pca.var,main="Scree plot", xlab="PC", ylab="Variation (%)")

pca.data = data.frame(Sample=rownames(pca$x),
                      X=pca$x[,1],
                      Y=pca$x[,2])
head(pca.data)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")
head(pca.data)
head(sub_data2)
pca.data.labeled = cbind(pca.data,sub_data2_unb[,'label'])
names(pca.data.labeled)[4] = "label"
head(pca.data.labeled)
ggplot(data=pca.data.labeled, aes(x=X, y=Y, colour=label))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")


ggplot(pca.data.labeled, aes(x=X, y=Y, colour=label)) +
  geom_point() + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_unb$label)
ggplot(pca.data.labeled, aes(x=X, y=Y)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_unb$label)

#n10
pca = prcomp(sub_data2_n10[,3:9], scale=TRUE)
plot(pca$x[,1], pca$x[,2])
pca.var = pca$sdev^2
pca.var.per = round(pca.var/sum(pca.var)*100,1)
barplot(pca.var,main="Scree plot", xlab="PC", ylab="Variation (%)")

pca.data = data.frame(Sample=rownames(pca$x),
                      X=pca$x[,1],
                      Y=pca$x[,2])
head(pca.data)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")
head(pca.data)
head(sub_data2)
pca.data.labeled = cbind(pca.data,sub_data2_n10[,'label'])
names(pca.data.labeled)[4] = "label"
head(pca.data.labeled)
ggplot(data=pca.data.labeled, aes(x=X, y=Y, colour=label))+
  geom_point()+#geom_text() +
  xlab(paste("PC1 - ",pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ",pca.var.per[2], "%", sep=""))+
  theme_bw()+
  ggtitle("PCA")


ggplot(pca.data.labeled, aes(x=X, y=Y, colour=label)) +
  geom_point() + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_n10$label)
ggplot(pca.data.labeled, aes(x=X, y=Y)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("PCA - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_n10$label)


###################################################
# t-SNE
head(iris)
iris_unique <- unique(iris) # Remove duplicates
iris_matrix <- as.matrix(iris_unique[,1:4])
set.seed(42) # Set a seed if you want reproducible results
tsne_out <- Rtsne(iris_matrix) # Run TSNE
plot(tsne_out$Y,col=iris_unique$Species)

head(sub_data2[,3:10])

mat = as.matrix(sub_data2[,3:9])
tsne_out <- Rtsne(mat, pca_scale=TRUE) # Run TSNE
plot(tsne_out$Y,col=sub_data2$label)
class(tsne_out$Y)
head(tsne_out$Y)
tsne.df=data.frame(tsne_out$Y)
tsne.df2 = cbind(tsne.df, sub_data2$label)
ggplot(tsne.df2, aes(x = X1, y = X2, colour=sub_data2$label)) +
  geom_point() + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2$label)
ggplot(tsne.df2, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2$label)

# unbalanced:
mat = as.matrix(sub_data2_unb[,3:9])
tsne_out <- Rtsne(mat, pca_scale=TRUE) # Run TSNE
plot(tsne_out$Y,col=sub_data2_unb$label+1)
class(tsne_out$Y)
head(tsne_out$Y)
tsne.df=data.frame(tsne_out$Y)
tsne.df2 = cbind(tsne.df, sub_data2_unb$label)
ggplot(tsne.df2, aes(x = X1, y = X2, colour=sub_data2_unb$label)) +
  geom_point() + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_unb$label)
ggplot(tsne.df2, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_unb$label)

#n10:
mat = as.matrix(sub_data2_n10[,3:9])
tsne_out <- Rtsne(mat, pca_scale=TRUE) # Run TSNE
plot(tsne_out$Y,col=sub_data2_n10$label)
class(tsne_out$Y)
head(tsne_out$Y)
tsne.df=data.frame(tsne_out$Y)
tsne.df2 = cbind(tsne.df, sub_data2_n10$label)
ggplot(tsne.df2, aes(x = X1, y = X2, colour=sub_data2_n10$label)) +
  geom_point() + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_n10$label)
ggplot(tsne.df2, aes(x = X1, y = X2)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("t-sne - f1chr4 - true and false edges") + xlab("x") + ylab("y")+ theme(text = element_text(size=16)) + facet_wrap( ~ sub_data2_n10$label)


###################################################
# plots for wrong edges:

ggplot(data, aes(x = n, y = w)) +
  geom_density_2d() + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("w")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(data, aes(x = n, y = w)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("w")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)


ggplot(data, aes(x = n, y = w_tfidf)) +
  geom_density_2d() + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("w_tfidf")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(data, aes(x = n, y = w_tfidf)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("w_tfidf")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)

# n vs n_jaccard
ggplot(data, aes(x = n, y = n_jaccard)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("n_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(data, aes(x = n, y = n_jaccard)) +
  geom_point(size=.1)+stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("n_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)

# n vs n_tfidf
ggplot(data, aes(x = n, y = n_tfidf)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("n_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(data, aes(x = n, y = n_tfidf)) +
  geom_point(size=.1)+stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("n") + ylab("n_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)

# w vs w_jaccard
ggplot(data, aes(x = w, y = w_jaccard)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("w") + ylab("w_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)
ggplot(data, aes(x = w, y = w_jaccard)) +
  geom_point(size=.1)+stat_density_2d(aes(fill = stat(level)), geom = "polygon") + ggtitle("f1chr4 - true and false edges stats") + xlab("w") + ylab("w_jaccard")+ theme(text = element_text(size=16)) + facet_wrap( ~ label)


#######################################################################
#### unlabeled file:



ad = "/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/f1chr4.k32-w32.n100-1000.c2-x.physlr.overlap.nn.edges.tsv"

data = read.table(ad, sep = "\t", header = TRUE)

# subsample 
#sub = data.frame(data_t[sample(nrow(data_t), dim(data_f)[1]), ])
#data = rbind(sub_t,data_f)
sub = data.frame(data[sample(nrow(data), 100000), ])

dim(data)
head(data)
tail(data)
class(data[1,3])
#data$label <- as.factor(data$label)

# hitogram of all features

# plot 2 features at a time
sub_data2 = data.frame(data[sample(nrow(data), 30000), ])
head(sub_data2)
ggplot(sub_data2, aes(x = n, y = sn1)) +
  geom_point() + ggtitle("f1chr4 n vs. ns1") + xlab("n") + ylab("ns1")+ theme(text = element_text(size=16))
ggplot(sub_data2, aes(x = n, y = sn1)) +
  geom_density_2d() + ggtitle("f1chr4 n vs. ns1") + xlab("n") + ylab("ns1")+ theme(text = element_text(size=16))
