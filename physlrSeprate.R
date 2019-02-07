library("ggplot2")
library("cowplot")
library("igraph")
library("lsa")


dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/AAGCCGCCAGGATCGA.tsv" # 38 HARD!!!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AATCGTGAGGAGTCTG.tsv" # HARD!!!

dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/AAACCTGCACGCTTTC.tsv" #22  #DONE poorly..!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/AACGGAGTCCTGCACT.tsv" #307 #DONE
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/AACTCCCTCTCCTGGT.tsv" #94 #DONE
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/AGGAGACTCGAAAGGC.tsv" #90 DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/ATTACTCAGTCGGGAT.tsv" #151 DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/CAAGGCCCATGAAGTA.tsv" #101 DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/GCGCGTAGTTTCCGGG.tsv" #103 DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/GGGACAAAGGACAGTC-GTACGTACAAAGGCTG.tsv" #84 #DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/GTTGCAACAACGTCTA.tsv" #184 #DONE
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/TCATTACGTCTCTGGG.tsv" #113 #DONE!
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/TCGGGCAAGAGCACCA.tsv" #118 #Done! #figures: dataSmilesOnUs_nodeDegreeHistogram
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/eg/TTTGTGTAGCCTGATT.tsv" #148 Done!

dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/GAGCTTACAGGATCTT.tsv" #simplDONE but No graphs available
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/AGTAGTCGTTCTCTCG-GTTTGTTCAGCGCATC.physlr.tsv" # No graphs available
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/CAAAGGGAGCACTCAT.tsv" # No files
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/CATATTCGTTATCACG.tsv" #No graphs available
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/GACTGCGCACCATGTA.tsv" #No graphs available
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/TCCGGTTAGTCCCGAC.tsv" #No graphs available

dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/TCGCTACAGGCCTCCA.tsv" #4?Weird file
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/CGAATGTTCCACCTTG.tsv" #2? file corrupted?
dat="/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/eg/CAGCTGGCACGGCTGT.tsv" #WHAT?

dat

#### Troubleshooting 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGTCATCGCTTAGA-1.tsv" #103
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGTCATCGCTTAGA-1_all.tsv" #532
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGTCACAGGTGTAG-1.tsv" # !!
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGTCACAGGTGTAG-1_all.tsv" #271
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AGCTGGCTCGCGCACA-1.tsv" # 91
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AGCTGGCTCGCGCACA-1_all.tsv" #615
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGCTATCGAGAGGT-1.tsv" # 57
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/TTTGCTATCGAGAGGT-1_all.tsv" #416
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTAGTAGGTCG-1.tsv" #191
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTAGTAGGTCG-1_all.tsv" #491

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTCAAGCGATG-1.tsv" #184 seems to be 1 dominant component
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTCAAGCGATG-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTGTATACACC-1.tsv" # 209
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTGTATACACC-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGCATCCACAGAT-1.tsv" #23 Seems to be 1 comp only
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGCATCCACAGAT-1_all.tsv" #304 strongly 1 comp
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTGTGGTCAGA-1.tsv" # too smaaaall
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTGTGGTCAGA-1_all.tsv" #206 strongly 1 comp
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGGTCACAACTTG-1.tsv" # 301 seem to be 1
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGGTCACAACTTG-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGGTTCTGGTGTA-1.tsv" #108
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGGTTCTGGTGTA-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGCAATGTAAG-1.tsv" # 97
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGCAATGTAAG-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGCACAGGAGT-1.tsv" #18
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGCACAGGAGT-1_all.tsv" #263 strongly 1 file

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGTCAAAGTAG-1.tsv" #170
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATGTGTCAAAGTAG-1_all.tsv" #
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACCAACAAGATCGGGT-1.tsv" #49 seem like a single comp
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACCAACAAGATCGGGT-1_all.tsv" #

######
#new lane of troubleshooting:

#sqCosine+:1 Tigmint:2 #Stronlgy seem to be 1 component!
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCTCTCTGCTG-1.tsv" # 40
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCTCTCTGCTG-1_all.tsv" #474 
#sqCosine+:1 Tigmint:2 #Stronlgy seem to be 1 component!
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCAGTGGACGT-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCAGTGGACGT-1_all.tsv" # 277

# sqCosine+:1 Tigmint:4
# NC_004354.4	15676986	15748255	AAACACCCAGCTATTG-1	204
# NT_037436.4	17640086	17674610	AAACACCCAGCTATTG-1	82
# NT_033777.3	30079806	30187837	AAACACCCAGCTATTG-1	294
# NC_004353.4	968505	1026671	AAACACCCAGCTATTG-1	144
#Stronlgy seem to be 1 component!
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCCAGCTATTG-1.tsv" # no such file
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/AAACACCCAGCTATTG-1_all.tsv" # 525

# sqCosine+:1 Tigmint:2
#NT_033779.5	13665321	13730955	ACATCTTCAAGCGATG-1	174
#NC_004353.4	176878	230807	ACATCTTCAAGCGATG-1	173
#Stronlgy seem to be 1 component!
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTCAAGCGATG-1.tsv" # 184
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/ACATCTTCAAGCGATG-1_all.tsv" # 493

###################################################################################
## new lane of comparison - after filtering to chromosome 4 and 99.54% consistency acheiving

##oversplitting:

#incon cuz of filtering resulting in a small poorly connected subgraph
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACCATCCAGGCTAGGT-1.tsv" # 9
#consistent! 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACCATCCAGGCTAGGT-1_all.tsv" # 608

#incon cuz of filtering resulting in a small poorly connected subgraph
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CTTAATCTCTTGTAGG-1.tsv" # 13
#consistent! 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CTTAATCTCTTGTAGG-1_all.tsv" # 918

#incon cuz of filtering resulting in a small poorly connected subgraph
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/TACTTCAAGGAATCCG-1.tsv" # 10
#consistent! 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/TACTTCAAGGAATCCG-1_all.tsv" # 402

#incon cuz of filtering resulting in a small poorly connected subgraph
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/TCACGAAAGTGTCCGC-1.tsv" # 25 
#consistent! 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/TCACGAAAGTGTCCGC-1_all.tsv" # 159

## undersplitting:

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/AATCCAGAGGTTGTAA-1.tsv" # 78
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/AATCCAGAGGTTGTAA-1_all.tsv" # 404

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACAGTGTAGCTAGTCT-1.tsv" # 87
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACAGTGTAGCTAGTCT-1_all.tsv" # 778

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACATCTTAGAGACTAT-1.tsv" # 151
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACATCTTAGAGACTAT-1_all.tsv" # 835

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACGTCAATCTGCGTAA-1.tsv" # 128
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACGTCAATCTGCGTAA-1_all.tsv" # 782

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACTGAGTAGCAGCCCT-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ACTGAGTAGCAGCCCT-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/AGTTGGTTCAAGCCGC-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/AGTTGGTTCAAGCCGC-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATAGTGCTCCTGTGGG-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATAGTGCTCCTGTGGG-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATATCGGTCGAATGTC-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATATCGGTCGAATGTC-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATATGCGAGACCCATT-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/ATATGCGAGACCCATT-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CAACTTTCAACCAGAG-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CAACTTTCAACCAGAG-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CACTTTATCAAACCGT-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CACTTTATCAAACCGT-1_all.tsv" # 

dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CGATCAAAGGTAGACC-1.tsv" # 
dat="/projects/btl/aafshinfard/projects/physlr-dev/data/subgraphs/weird/CGATCAAAGGTAGACC-1_all.tsv" # 

1 CGGATTAGTGCCACCA-1
1 CGTTCCATCTCAAGTG-1
1 CTACACCTCGCTGATA-1
1 CTGATCCGTAGGACTG-1
1 CTGTGCTCAGTGCGTC-1
1 GAGTGAGAGAGACTTA-1
1 GCAACCGGTTTAGGAA-1
1 GCGTGTGCAAACAACA-1
1 GCTATAGTCTATCCTA-1
1 GGCTTAAGTGTGCCTG-1
1 GTACGGCTCTCAAGTG-1
1 GTTGCAACAACGTCTA-1
1 TAAGGCTCACCAGTGC-1
1 TAGAAGAGTTGACATC-1
1 TAGACACCAACCAGAG-1
1 TCAACGACAAGGACTG-1
1 TTAGGACTCAGTGTGT-1
1 TTCCCAGGTGTCCCTT-1
1 TTCGAAGCACGTCGGT-1
1 TTCTCAACAGAGCCAA-1
1 TTGTTTGTCCCAAGAT-1

########################################
### Summary
a = read.table(dat, header = FALSE, col.names = paste0("V",seq_len(3)),as.is = "V3", fill = TRUE)
a = a[118:dim(a)[1],]
a[781,]
a[,3] = as.numeric(a[,3])
gr <- graph.data.frame(a)
gr = as.undirected(gr)
adj = as_adjacency_matrix(gr, type = c("both"), attr = NULL, edges = FALSE, names = FALSE, sparse=FALSE)
# Remove nodes with no edge to source node
dim(adj)
head(adj)
adj_filt = adj
sort(colSums(adj_filt))
length(colSums(adj_filt))

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
length(colSums(adj_filt))
#hist(colSums(adj_filt))
ggplot(as.data.frame(colSums(adj_filt)), aes(x=colSums(adj_filt)) ) + geom_histogram( binwidth= 3, position="identity" )
adj_orig = adj_filt
roll <- function( x , n ){
  if( n == 0 )
    return( x )
  c( tail(x,n) , head(x,-n) )
}
#which(colSums(adj_filt) == 28)
#filt_index = c(1:7,9:13,15:34)
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

#filt_index = c(1:(which(colSums(adj_filt)==28)-1),(which(colSums(adj_filt)==28)+1):34) 
#adj_filt = adj_filt[filt_index,filt_index]
#filt_index = c(1:(which(colSums(adj_filt)==25)-1),(which(colSums(adj_filt)==25)+1):33) 
#adj_filt = adj_filt[filt_index,filt_index]
#sort(colSums(adj_filt[filt_index,filt_index]))
#length(colSums(adj_filt[filt_index,filt_index]))

sort(colSums(adj_filt))
length(colSums(adj_filt))
#hist(colSums(adj_filt))
ggplot(as.data.frame(colSums(adj_filt)), aes(x=colSums(adj_filt)) ) + geom_histogram(  binwidth= 4, position="identity" )

adj_sq = adj_filt%*%adj_filt
ggplot(as.data.frame(as.vector(adj_sq)), aes(x=as.vector(adj_sq)) ) + geom_histogram( aes(y=..density..), position="identity" )
ggplot(as.data.frame(colSums(adj_sq)), aes(x=colSums(adj_sq)) ) + geom_histogram( aes(y=..density..), position="identity" )

cos = cosine(adj_filt)
cos[cos == "NaN"] = 0
#hist(cos)
ggplot(as.data.frame(as.vector(cos)), aes(x=as.vector(cos)) ) + geom_histogram( aes(y=..density..), binwidth= 0.05, position="identity" )
cos2 = cosine(adj_sq)
cos3 = cos2
cos3[cos3 == "NaN"] = 0
#hist((cos3),breaks = 40)
ggplot(as.data.frame(as.vector(cos3)), aes(x=as.vector(cos3)) ) + geom_histogram( aes(y=..density..), binwidth= 0.02, position="identity" )

dim(adj_orig)
plot(graph_from_adjacency_matrix(adj_orig))
plot(graph_from_adjacency_matrix(adj_filt))
plot(graph_from_adjacency_matrix(adj_sq))
adj_t = adj_filt
adj_t[cos[,] < 0.84] = 0
adj_t[cos3[,] < 0.88] = 0
sum(cos3[,] < 0.9)/sum((cos3[,] > -1))
sum(adj_filt != 0)
sum(adj_t != adj_filt)
grnew = graph_from_adjacency_matrix(adj_t)
comps = components(grnew)
plot(grnew)

### using cos
remStat = as.data.frame({})
j = 1
adj_t = adj_filt
for(i in 0:100){
  print(i)
  adj_t[cos[,] < i/100] = 0
  remStat[j,1] = sum(adj_t != adj_filt)
  remStat[j,2] = i/100
  grnew = graph_from_adjacency_matrix(adj_t)
  comps = components(grnew)
  remStat[j,3] = sum(comps$csize > 1)
  remStat[j,4] = sum(comps$csize > 0)
  if( j>1 ){
    if(remStat[j,4] == remStat[j-1,4] & remStat[j,3] == remStat[j-1,3])
      remStat[j,5] = ""
    else
      remStat[j,5] = paste(as.character(remStat[j,3]),"-",as.character(remStat[j,4]),sep = "")
  }else
    remStat[j,5] = paste(as.character(remStat[j,3]),"-",as.character(remStat[j,4]),sep = "")
  j=j+1
}

names(remStat) <-c("V1","V2","ConComps","ConComps2","mixLabel")
head(remStat)
p1 <- ggplot(remStat, aes(x=V2,y=V1,label=mixLabel) ) + geom_col() + 
  geom_text(aes(label=mixLabel),hjust=0, col="red", vjust=0.5, angle = 90) +
  ggtitle("Threshold's effect on edge removal") +
  labs(x= "Threshold",y= "Number of edges removed") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1)) + xlim(0.4,1.01) + ylim(0,max(remStat[,1])+(max(remStat[,1])/10) )

p1


### Using cos3
remStat3 = as.data.frame({})
j = 1
adj_t = adj_filt
for(i in 0:100){
  print(i)
  adj_t[cos3[,] < i/100] = 0
  remStat3[j,1] = sum(adj_t != adj_filt)
  remStat3[j,2] = i/100
  grnew = graph_from_adjacency_matrix(adj_t)
  comps = components(grnew)
  remStat3[j,3] = sum(comps$csize > 1)
  remStat3[j,4] = sum(comps$csize > 0)
  if( j>1 ){
    if(remStat3[j,4] == remStat3[j-1,4] & remStat3[j,3] == remStat3[j-1,3])
    remStat3[j,5] = ""
    else
      remStat3[j,5] = paste(as.character(remStat3[j,3]),"-",as.character(remStat3[j,4]),sep = "")
  }else
    remStat3[j,5] = paste(as.character(remStat3[j,3]),"-",as.character(remStat3[j,4]),sep = "")
  j=j+1
}

names(remStat3) <-c("V1","V2","ConComps","ConComps2","mixLabel")
head(remStat3)
p1 <- ggplot(remStat3, aes(x=V2,y=V1,label=mixLabel) ) + geom_col() + 
  geom_text(aes(label=mixLabel),hjust=0, col="red", vjust=0.5, angle = 90) +
  ggtitle("Threshold's effect on edge removal") +
  labs(x= "Threshold",y= "Number of edges removed") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1)) + xlim(0,1.01) + ylim(0,max(remStat3[,1])+(max(remStat3[,1])/10) )
p1
#p1 <- ggplot(remStat, aes(x=V2,y=V1) ) + geom_col() + ggtitle("Threshold's effect on edge removal") +
#  labs(x= "Threshold",y= "Number of edges removed") + 
#  theme(text = element_text(size=14),
#       axis.text.x = element_text(angle=0, hjust=1)) + xlim(0.5,1.01)
  #ylim( 0,rem[(dim(rem)[1]-1),1]*2 )
#remStatM = melt(remStat[,2:4], id.vars="V2")
#p2 <- ggplot(remStatM, aes(V2,value, col=variable)) + geom_col(position = "dodge") + ggtitle("Threshold's effect on edge removal") +
#  labs(x= "Threshold",y= "Number of edges removed") + 
#  theme(text = element_text(size=14),
#        axis.text.x = element_text(angle=0, hjust=1)) + xlim(0.5,1.01)
#legend <- get_legend(p2)
#grid.arrange(p1, p2, legend, ncol = 1, heights = c(1, 1))
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
adj_t[cos3[,] < 0.97 ] = 0
sum(adj_t != adj_filt | adj_t == adj_filt)
sum(adj_t != adj_filt)
grnew = graph_from_adjacency_matrix(adj_t)
plot(grnew)

comps = components(grnew)
comps


##### extra codes:
library(FinCal)
library(reshape2)

dat <- get.ohlc.yahoo('AAPL', '2015-12-01', '2015-12-31')
#dat$date <- strptime(dat$date, "%Y-%m-%d")
dat$date <- as.Date(dat$date, "%Y-%m-%d")
dat$times <- seq(nrow(dat))
mm <- melt(subset(dat, select=c(times,adjusted, volume)), id.var="times")

