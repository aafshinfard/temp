dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AAACCTGCACGCTTTC.tsv" #22 #DONE poorly..!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/TCATTACGTCTCTGGG.tsv" #113 #DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/GGGACAAAGGACAGTC-GTACGTACAAAGGCTG.tsv" #84 #DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AACTCCCTCTCCTGGT.tsv" #94 #DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AAGCCGCCAGGATCGA.tsv" ## HARD!!!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AATCGTGAGGAGTCTG.tsv" ## HARD!!!

########################################
### Summary
a = read.table(dat, header = FALSE, col.names = paste0("V",seq_len(3)),as.is = "V3", fill = TRUE)
a = a[22:dim(a)[1],]
a[,3] = as.numeric(a[,3])
gr <- graph.data.frame(a)
gr = as.undirected(gr)
adj = as_adjacency_matrix(gr, type = c("both"), attr = NULL, edges = FALSE, names = FALSE, sparse=FALSE)
# remove source node  
adj_filt = adj
if(max(colSums(adj) )== dim(adj)[1]-1 ){
  filt_index ={}
  if(which.max(colSums(adj))>1)
    filt_index = c(filt_index,1:( which.max(colSums(adj)) -1))
  if(which.max(colSums(adj))<dim(adj)[1])
    filt_index = c(filt_index,( which.max(colSums(adj)) +1):dim(adj)[1])
  adj_filt = adj[filt_index,filt_index]
}
# remove highly connected 'outlier' node
roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}
rolled = roll(sort(colSums(adj_filt)),1) - sort(colSums(adj_filt))
rolled = rolled[2:length(colSums(adj_filt))]
if(min(rolled)==rolled[length(rolled)]){
  filt_index ={}
  if(which.max(colSums(adj_filt))>1)
    filt_index = c(filt_index,1:( which.max(colSums(adj_filt)) -1))
  if(which.max(colSums(adj_filt))<dim(adj_filt)[1])
    filt_index = c(filt_index,( which.max(colSums(adj_filt)) +1):dim(adj_filt)[1])
  adj_filt = adj_filt[filt_index,filt_index]
}
adj_sq = adj_filt%*%adj_filt

cos = cosine(adj_filt)
cos[cos == "NaN"] = 0
hist(cos)
cos2 = cosine(adj_sq)
cos3 = cos2
cos3[cos3 == "NaN"] = 0
hist(cos3)

plot(graph_from_adjacency_matrix(adj_filt))
plot(graph_from_adjacency_matrix(adj_sq))
adj_t = adj_filt
adj_t[cos[,] < 0.7] = 0
#adj_t[cos3[,] < 0.98 ] = 0
sum(adj_t != adj_filt | adj_t == adj_filt)
sum(adj_t != adj_filt)
grnew = graph_from_adjacency_matrix(adj_t)
plot(grnew)

########################################
##### FULL

a = read.table(dat, header = FALSE, col.names = paste0("V",seq_len(3)),as.is = "V3", fill = TRUE)
#a = read.graph(format = "dl", "/projects/btl_scratch/sjackman/physlr/eg/TCATTACGTCTCTGGG.gv")
#a = read.graph( "/projects/btl_scratch/sjackman/physlr/eg/TCATTACGTCTCTGGG.tsv")
a = a[94:dim(a)[1],]
a[,3] = as.numeric(a[,3])
gr <- graph.data.frame(a)
#adjdir = as_adjacency_matrix(gr, type = c("both"), attr = NULL, edges = FALSE, names = FALSE, sparse=FALSE)
gr = as.undirected(gr)
adj = as_adjacency_matrix(gr, type = c("both"), attr = NULL,
                          edges = FALSE, names = FALSE, sparse=FALSE)
hist(colSums(adj))
sort(colSums(adj))
(dim(adj))
plot(gr)
#OutVals = boxplot(colSums(adj))$out
#which(x %in% OutVals)
#remove_outliers <- function(x, na.rm = TRUE, ...) {
#  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
#  H <- 1.5 * IQR(x, na.rm = na.rm)
#  y <- x
#  y[x < (qnt[1] - H)] <- NA
#  y[x > (qnt[2] + H)] <- NA
#  y
#}
# remove source node  
adj_filt = adj
if(max(colSums(adj) )== dim(adj)[1]-1 ){
  filt_index ={}
  if(which.max(colSums(adj))>1)
    filt_index = c(filt_index,1:( which.max(colSums(adj)) -1))
  if(which.max(colSums(adj))<dim(adj)[1])
    filt_index = c(filt_index,( which.max(colSums(adj)) +1):dim(adj)[1])
  adj_filt = adj[filt_index,filt_index]
}
# remove highly connected 'outlier' node
sort(colSums(adj_filt))
roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}
rolled = roll(sort(colSums(adj_filt)),1) - sort(colSums(adj_filt))
rolled = rolled[2:length(colSums(adj_filt))]
if(min(rolled)==rolled[length(rolled)]){
  filt_index ={}
  if(which.max(colSums(adj_filt))>1)
    filt_index = c(filt_index,1:( which.max(colSums(adj_filt)) -1))
  if(which.max(colSums(adj_filt))<dim(adj_filt)[1])
    filt_index = c(filt_index,( which.max(colSums(adj_filt)) +1):dim(adj_filt)[1])
  adj_filt = adj_filt[filt_index,filt_index]
}
sort(colSums(adj_filt))

hist(colSums(adj_filt))
sort(colSums(adj_filt))
(dim(adj_filt))
adj_sq = adj_filt%*%adj_filt
hist(colSums(adj_sq[,]>0))
sum(adj_sq == 0)
sort(colSums(adj_sq))

#dim(adj_filt)
#dim(adj_sq)
cos = cosine(adj_filt)
cos[cos == "NaN"] = 0
#cosdir = cosine(adjdir)
#cosdir[cosdir == "NaN"] = 0
hist(cos)
#hist(cosdir)
cos2 = cosine(adj_sq)
cos3 = cos2
cos3[cos3 == "NaN"] = 0
hist(cos3)
sum(adj == 0 )/ sum(adj == 0 | adj > 0)
sum(cos3 == 0 )/ sum(cos3 == 0 | cos3 > 0)

adj_sq = adj_filt%*%adj_filt
plot(graph_from_adjacency_matrix(adj_filt))
plot(graph_from_adjacency_matrix(adj_sq))
adj_t = adj_filt
#adj_t[cos[,] < 0.7] = 0
adj_t[cos3[,] < 0.7 ] = 0
sum(adj_t != adj_filt | adj_t == adj_filt)
sum(adj_t != adj_filt)
grnew = graph_from_adjacency_matrix(adj_t)
plot(grnew)

comps = components(grnew)
comps
