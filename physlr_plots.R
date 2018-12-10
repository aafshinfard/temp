library(ggplot2)

##############################################################################
############# Physlr node edge data 

addnodeg = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/run1/f1chr4.physlr.overlap.neighbor.stat.tsv"
addnodeg = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/prevRuns/run5/f1chr4.physlr.overlap.n50.mol.tsv"
addnodeg = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr/f1chr4.physlr.overlap.neighbor.stat.tsv"
addnodeg = "/projects/btl/aafshinfard/projects/physlr/physlr-cloned/physlr-old/prevRuns/simMat/f1chr4.physlr.overlap.neighbor.noFilter_simMat_Histo.tsv"

### If it's Histogram of cosine sim:
# go to the end of this code


### If it's |E| vs |V| :
nodeg = 
  read.table(addnodeg, sep = "\t", header = FALSE)
nodeg
head(nodeg)
tail(nodeg)
sum(nodeg$V1==0)
sum(nodeg$V1==1)
sum(nodeg$V1==1)/sum(nodeg$V1==0)


sum(nodeg$V3==0)
nodeg = nodeg[nodeg$V3!=0,]
nodeg = nodeg[nodeg$V2!=0,]
sum(nodeg$V2==0)

#plot1 = ggplot(data = hist10X, aes(x = V1, y=V2)) + geom_jitter(shape=21, width=1, height=.5, fill="#A1C0E5", color="#639DD2") + 
 # theme_minimal()
#plot1
#plot1 = ggplot(data = nodeg, aes(x = V1, y=V2)) + geom_jitter(shape=21, width=1, height=.5, fill="#A1C0E5") + 
#  theme_minimal()
#plot1

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
vv_plot+geom_step(directon="hv",aes(x=1:dim(vv)[1], y=as.vector(t(vv))))+xlim(c(0,120))+ylim(c(0,300000000))
