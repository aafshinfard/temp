library(ggplot2)

##############################################################################
############# Physlr node edge data 

#### Fly:
add_nodeg0 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/f1.n100-2000.physlr.overlap.subgraphs_stats"
title="Fly - whole data - n0 - subsampled to 30K subgraphs"
add_nodeg20 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n20.subgraphs_stats"
title="Fly - whole data - n20 - subsampled to 30K subgraphs"
add_nodeg50 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n50.subgraphs_stats"
title="Fly - whole data - n50 - subsampled to 30K subgraphs"
add_nodeg70 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n70.subgraphs_stats"
title="Fly - whole data - n70 - subsampled to 30K subgraphs"
add_nodeg100 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n100.subgraphs_stats"
title="Fly - whole data - n100 - subsampled to 30K subgraphs"
add_nodeg200 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n200.subgraphs_stats"
title="Fly - whole data - n200 - subsampled to 30K subgraphs"
### If it's Histogram of cosine sim:
# go to the end of this code

### Fish:
add_nodeg100 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n100.subgraphs_stats"


add_nodeg20 = "/projects/btl_scratch/aafshinfard/phys-dev/physlr/data/fish.indexlr.n100-2000.c2-x.physlr.overlap.n20.subgraphs_stats"
add_nodeg20_mol = "/projects/btl_scratch/aafshinfard/phys-dev/physlr/data/fish.indexlr.n100-2000.c2-x.physlr.overlap.n20.mol.subgraphs_stats"

### If it's |E| vs |V| :
nodeg0 = 
  read.table(add_nodeg0, sep = "\t", header = FALSE, skip = 1)
nodeg20 = 
  read.table(add_nodeg20, sep = "\t", header = FALSE, skip = 1)
nodeg20_mol = 
  read.table(add_nodeg20_mol, sep = "\t", header = FALSE, skip = 1)
nodeg50 = 
  read.table(add_nodeg50, sep = "\t", header = FALSE, skip = 1)
nodeg70 = 
  read.table(add_nodeg70, sep = "\t", header = FALSE, skip = 1)
nodeg100 = 
  read.table(add_nodeg100, sep = "\t", header = FALSE, skip = 1)
nodeg200 = 
  read.table(add_nodeg200, sep = "\t", header = FALSE, skip = 1)

#nodeg
names(nodeg0)=c("Barcodes","nodes","edges","alpha")
names(nodeg20)=c("Barcodes","nodes","edges","alpha")
names(nodeg20_mol)=c("Barcodes","nodes","edges","alpha")
names(nodeg50)=c("Barcodes","nodes","edges","alpha")
names(nodeg70)=c("Barcodes","nodes","edges","alpha")
names(nodeg100)=c("Barcodes","nodes","edges","alpha")
names(nodeg200)=c("Barcodes","nodes","edges","alpha")
head(nodeg0)
tail(nodeg0)

#head(nodeg[,"nodes","edges"])
#headi = head(nodeg)
#headi[,2:3]

a=nodeg_filt0
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt0 = nodeg0[nodeg0[,"nodes"]>th,]

a=nodeg20
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt20 = nodeg20[nodeg20[,"nodes"]>th,]

a=nodeg20_mol
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt20_mol = nodeg20_mol[nodeg20_mol[,"nodes"]>th,]


a=nodeg_filt50
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt50 = nodeg50[nodeg50[,"nodes"]>th,]

a=nodeg_filt70
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt70 = nodeg70[nodeg70[,"nodes"]>th,]

a=nodeg_filt100
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt100 = nodeg100[nodeg100[,"nodes"]>th,]

a=nodeg_filt200
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt200 = nodeg200[nodeg200[,"nodes"]>th,]

#dim(nodeg)[1]
dim(nodeg_filt0)[1]

#head(nodeg_filt)
maxi = 70000
nodeg_filt2_0 = nodeg_filt0[1:maxi,]
nodeg_filt2_20 = nodeg_filt20[1:maxi,]
nodeg_filt2_20_mol = nodeg_filt20_mol[1:maxi,]
nodeg_filt2_50 = nodeg_filt50[1:maxi,]
nodeg_filt2_70 = nodeg_filt70[1:maxi,]
nodeg_filt2_100 = nodeg_filt100[1:maxi,]
nodeg_filt2_200 = nodeg_filt200[1:maxi,]

title="Fish - after mol2bar - n20 - subsampled to 70K subgraphs"
#ploti = ggplot(nodeg, aes(x=`nodes`, y=`edges`,color=`multip-comp`)) + geom_point()+
ploti = ggplot(nodeg_filt2_20_mol, aes(x=`nodes`, y=`edges`)) + geom_point(size=0.001)+
  ggtitle(title) +
  labs(x= "Number of nodes",y= "Number of edges") + 
  xlim(0,300)+
  ylim(0,8000)+#ylim(0,2.5*10^5)+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))
ploti
ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`alpha`)) + geom_point(size=0.001)+
  ggtitle(title) +
  labs(x= "Number of nodes",y= "Number of edges") + 
  xlim(0,300)+
  #ylim(0,8000)+#ylim(0,2.5*10^5)+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))+
  geom_point(data = nodeg_filt2_20_mol,size=0.001,color="red")
ploti


1
2
3
4
5
6
7
8
ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`alpha`)) +
  geom_hex(bins=100) +
  theme_bw()
ploti

ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`alpha`)) +
  geom_density_2d() +
  theme_bw()
ploti

ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`alpha`)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
  scale_fill_viridis_c()+
  theme_bw()
ploti

ploti =ggplot(nodeg_filt0, aes(x=`nodes`, y=`edges`)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )
ploti
####################################################
########################## Tigmint coloring
### Fly:
add_nodeg0 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/f1.n100-2000.physlr.overlap.subgraphs_stats.molCount"
title="Fly - whole data - n0 - subsampled to 30K subgraphs"
add_nodeg20 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n20.subgraphs_stats.molCount"
title="Fly - whole data - n20 - subsampled to 30K subgraphs"
add_nodeg50 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n50.subgraphs_stats.molCount"
title="Fly - whole data - n50 - subsampled to 30K subgraphs"
add_nodeg70 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n70.subgraphs_stats."
title="Fly - whole data - n70 - subsampled to 30K subgraphs"
add_nodeg100 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n100.subgraphs_stats.molCount"
title="Fly - whole data - n100 - subsampled to 30K subgraphs"
add_nodeg200 = "/projects/btl/aafshinfard/projects/physlr-dev/data/f1.n100-2000.physlr.overlap.n200.subgraphs_stats.molCount"
title="Fly - whole data - n200 - subsampled to 30K subgraphs"
### If it's Histogram of cosine sim:
# go to the end of this code
### Fish:
add_nodeg0 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.subgraphs_stats.molCount"
add_nodeg20 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n20.subgraphs_stats.molCount"
add_nodeg29 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n29.subgraphs_stats.molCount"
add_nodeg50 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n40.subgraphs_stats.molCount"
add_nodeg53 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n53.subgraphs_stats.molCount"
add_nodeg100 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n100.subgraphs_stats.molCount"
add_nodeg200 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n200.subgraphs_stats.molCount"

add_nodeg0 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.subgraphs_stats.molCount"
add_nodeg23 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n23.subgraphs_stats.molCount"
add_nodeg41 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n41.subgraphs_stats.molCount"
add_nodeg65 = "/projects/btl_scratch/aafshinfard/physlr-dev/data/fish.indexlr.n100-5000.c2-x.physlr.overlap.n65.subgraphs_stats.molCount"

### 
nodeg0 = 
  read.table(add_nodeg0, sep = " ", header = FALSE, skip = 1)
nodeg20 = 
  read.table(add_nodeg20, sep = " ", header = FALSE, skip = 1)
nodeg29 = 
  read.table(add_nodeg29, sep = " ", header = FALSE, skip = 1)
nodeg50 = 
  read.table(add_nodeg50, sep = " ", header = FALSE, skip = 1)
nodeg53 = 
  read.table(add_nodeg53, sep = " ", header = FALSE, skip = 1)
nodeg70 = 
  read.table(add_nodeg70, sep = " ", header = FALSE, skip = 1)
nodeg100 = 
  read.table(add_nodeg100, sep = " ", header = FALSE, skip = 1)
nodeg200 = 
  read.table(add_nodeg200, sep = " ", header = FALSE, skip = 1)



nodeg0 = 
  read.table(add_nodeg0, sep = " ", header = FALSE, skip = 1)
nodeg23 = 
  read.table(add_nodeg23, sep = " ", header = FALSE, skip = 1)
nodeg41 = 
  read.table(add_nodeg41, sep = " ", header = FALSE, skip = 1)
nodeg65 = 
  read.table(add_nodeg65, sep = " ", header = FALSE, skip = 1)
  
#nodeg
names(nodeg0)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg20)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg29)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg50)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg53)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg70)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg100)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg200)=c("Barcodes","nodes","edges","alpha","molCount")
head(nodeg100)
tail(nodeg100)


names(nodeg0)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg23)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg41)=c("Barcodes","nodes","edges","alpha","molCount")
names(nodeg65)=c("Barcodes","nodes","edges","alpha","molCount")

#head(nodeg[,"nodes","edges"])
#headi = head(nodeg)
#headi[,2:3]

a=nodeg0
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt0 = nodeg0[nodeg0[,"nodes"]>th,]

a=nodeg_filt20
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt20 = nodeg20[nodeg20[,"nodes"]>th,]

a=nodeg29
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt29 = nodeg29[nodeg29[,"nodes"]>th,]

a=nodeg_filt50
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt50 = nodeg50[nodeg50[,"nodes"]>th,]

a=nodeg_filt53
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt53 = nodeg53[nodeg53[,"nodes"]>th,]

a=nodeg_filt70
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt70 = nodeg70[nodeg70[,"nodes"]>th,]

a=nodeg_filt100
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt100 = nodeg100[nodeg100[,"nodes"]>th,]

a=nodeg_filt200
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt200 = nodeg200[nodeg200[,"nodes"]>th,]



a=nodeg0
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt0 = nodeg0[nodeg0[,"nodes"]>th,]

a=nodeg23
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt23 = nodeg23[nodeg23[,"nodes"]>th,]

a=nodeg41
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt41 = nodeg41[nodeg41[,"nodes"]>th,]

a=nodeg65
th=sort((a[,"nodes"]))[length(a[,"nodes"])/10]
nodeg_filt65 = nodeg65[nodeg65[,"nodes"]>th,]

#dim(nodeg)[1]
dim(nodeg_filt0)[1]

#head(nodeg_filt)
maxi = 70000
nodeg_filt2_0 = nodeg_filt0[1:maxi,]
nodeg_filt2_20 = nodeg_filt20[1:maxi,]
nodeg_filt2_29 = nodeg_filt29[1:maxi,]
nodeg_filt2_50 = nodeg_filt50[1:maxi,]
nodeg_filt2_53 = nodeg_filt53[1:maxi,]
nodeg_filt2_70 = nodeg_filt70[1:maxi,]
nodeg_filt2_100 = nodeg_filt100[1:maxi,]
nodeg_filt2_200 = nodeg_filt200[1:maxi,]

nodeg_filt2_0 = nodeg_filt0[1:maxi,]
nodeg_filt2_23 = nodeg_filt23[1:maxi,]
nodeg_filt2_41 = nodeg_filt41[1:maxi,]
nodeg_filt2_65 = nodeg_filt65[1:maxi,]





nodeg_filt2_0 = nodeg_filt2_0[nodeg_filt2_0[,"molCount"]<9,]
nodeg_filt2_20 = nodeg_filt2_20[nodeg_filt2_20[,"molCount"]<9,]
nodeg_filt2_29 = nodeg_filt2_29[nodeg_filt2_29[,"molCount"]<9,]
nodeg_filt2_50 = nodeg_filt2_50[nodeg_filt2_50[,"molCount"]<9,]
nodeg_filt2_53 = nodeg_filt2_53[nodeg_filt2_53[,"molCount"]<9,]
nodeg_filt2_100 = nodeg_filt2_100[nodeg_filt2_100[,"molCount"]<9,]
nodeg_filt2_200 = nodeg_filt2_200[nodeg_filt2_200[,"molCount"]<7,]

nodeg_filt2_0 = nodeg_filt2_0[nodeg_filt2_0[,"molCount"]<9,]
nodeg_filt2_23 = nodeg_filt2_23[nodeg_filt2_23[,"molCount"]<9,]
nodeg_filt2_41 = nodeg_filt2_41[nodeg_filt2_41[,"molCount"]<9,]
nodeg_filt2_65 = nodeg_filt2_65[nodeg_filt2_65[,"molCount"]<9,]

#ploti = ggplot(nodeg_filt2_0, aes(x=`nodes`, y=`edges`)) +
#  geom_bin2d(bins = 400) +
#  theme_bw()

#ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`edges`)) +
#  stat_density_2d(aes(fill = ..level..), geom = "polygon")
#ploti

bb=nodeg_filt2_20
head(bb)

title="Fish - whole data - n65 - subsampled to 70K subgraphs"
#ploti = ggplot(nodeg_filt2_20, aes(x=`nodes`, y=`edges`,color=nodeg_filt2_20[,"molCount"])) + geom_point(size=0.001)+
ploti = ggplot(nodeg_filt2_65, aes(x=`nodes`, y=`edges`)) + geom_point(size=0.001,colour=nodeg_filt2_65[,"molCount"])+
  ggtitle(title) +
  labs(x= "Number of nodes",y= "Number of edges") + 
  #xlim(0,1300)+
  #ylim(0,2.5*10^5)+
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))+
  scale_color_gradientn(colours = rainbow(5))
ploti


a
a
a

head(nodeg_filt0[,])
nodeg0_1d = nodeg_filt0[nodeg_filt0[,"nodes"]<810,]
nodeg0_1d = nodeg0_1d[nodeg0_1d[,"nodes"]>800,]
head(nodeg0_1d)
hist(nodeg0_1d[,"edges"],breaks=30)
ggplot(as.data.frame(as.vector(nodeg0_1d[,"edges"])), aes(x=as.vector(nodeg0_1d[,"edges"])) ) + geom_histogram( aes(y=..density..), binwidth= 0.05, position="identity" )
plot1


plot1 = ggplot(nodeg, aes(x=`|V|`, y=`|E|`)) + geom_point(color="#639DD2")+
  ggtitle("Fly chromosome 4") +
  labs(x= "Number of vertices |V|",y= "Number of edges |E|") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))

plot1


##############################################################################
##############################################################################
################################## OLD CODES:
#nodeg0 = nodeg[nodeg['V1']==0,]
#nodeg1 = nodeg[nodeg['V1']==1,]

names(nodeg)=c("multip-comp","|V|","|E|")
nodeg$`multip-comp` <- as.factor(nodeg$`multip-comp`)

head(nodeg)
ploti = ggplot(nodeg, aes(x=`|V|`, y=`|E|`,color=`multip-comp`)) + geom_point()+
  ggtitle("Fly chromosome 4") +
  labs(x= "Number of vertices |V|",y= "Number of edges |E|") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))

ploti
plot1 = ggplot(nodeg, aes(x=`|V|`, y=`|E|`)) + geom_point(color="#639DD2")+
  ggtitle("Fly chromosome 4") +
  labs(x= "Number of vertices |V|",y= "Number of edges |E|") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1))

plot1
egCount <- function(x,a){a*(((x^2)/2)-(x/2))}
egCount2 <- function(x,a){a*(((x^2)/4)-(x/2))}
#plot(nodeg$`|V|`,nodeg$`|E|`)
#curve(egCount(x,1),add=TRUE)
c = nodeg$`|E|`
d = nodeg$`|V|`
Zfit1 <- nls(c~egCount(d,a),start=list(a=0.5))
Zfit1

nodeg1c = nodeg[nodeg[,1]==0,]
c = nodeg1c$`|E|`
d = nodeg1c$`|V|`
Zfit2 <- nls(c~egCount(d,a),start=list(a=0.5))
Zfit2

nodeg2c = nodeg[nodeg[,1]==1,]
c = nodeg2c$`|E|`
d = nodeg2c$`|V|`
Zfit3 <- nls(c~egCount(d,a),start=list(a=0.5))
Zfit3

#curve(egCount(x,0.3916),add=TRUE)

#dummyknot = rep(0,length(nodeg$`|V|`))
#dummyknot[nodeg$`|V|`>25]=1
#nodeg$`|V|25` = nodeg$`|V|2` - 25

#fun.n <- function(x) x
alpha = 0.45
e_alpha = 0.5978
fun.1 <- function(x) alpha*((x^2-x)/2)
fun.2 <- function(x) alpha*((x^2-(2*x))/4)
fun.estimated <- function(x) e_alpha*((x^2-x)/2)
e_alpha2 = 0.5785
fun.estimated2 <- function(x) e_alpha2*((x^2-x)/2)

ploti + stat_function(fun = fun.estimated,aes(color="fun")) + scale_color_manual(name = "fun", values = c("red","blue","red"))
ploti + layer(geom = "path",        # Default. Can be omitted.
              stat = "function",
              fun = fun.estimated,          # Give function
              mapping = aes(color = "fun.estimated") # Give a meaningful name to color
) +
  scale_color_manual(name = "Function", values = c("red"))


##############################################################################
############# Physlr Similarity Matrix


add = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/simMat/f1chr4.physlr.overlap.neighbor_simMat_Histo.tsv"
add = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/simMat/f1chr4.physlr.overlap.neighbor.n50_simMat_Histo_diagSet1.tsv"
add = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/simMat/f1chr4.physlr.overlap.neighbor.n50_simMat_Histo.tsv"
add = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/simMat/f1chr4.physlr.overlap.neighbor.n20_simMat_Histo.tsv"
add = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/prevRuns/simMat/f1chr4.physlr.overlap.neighbor.squared_Histo.tsv"

vv = read.table(add, sep = "\t", header = FALSE)
vv_plot = ggplot(vv)
vv_plot+geom_step(directon="hv",aes(x=1:dim(vv)[1], y=as.vector(t(vv))))+xlim(c(0,120))+ylim(c(0,120000000))+
  ggtitle("Threshold's effect on edge removal - over all subgraphs") +
  labs(x= "Threshold X 100",y= "Number of edges removed") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(angle=0, hjust=1)) 
