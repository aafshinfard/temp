# Thank Lauren for the initial version of the script

library(ggplot2)
library(dplyr)

mol <- read.csv("/projects/btl/lcoombe/git/physlr/data/zebrafish/zebrafish.zebrafish_lr.a0.65.d50000.n5.q1.s2000.molecule.bed",
                sep="\t", header=F)
colnames(mol) <- c("ref", "start", "end", "bx", "num_reads")

mol <- mol %>% mutate(Coverage=(num_reads*150)/(end-start))
mol <- mol %>% mutate(Length=(end-start))
mol2 = mol[1:100000,]
ggplot(mol, aes(x=Coverage)) +
  geom_histogram(stat="bin", binwidth=0.001)+
  scale_x_continuous(limits=c(0, 1))

ggplot(mol2, aes(x=Length)) +
  geom_histogram(stat="bin", binwidth=1000)+xlim(0,200000)

ggplot(mol2, aes(x=num_reads)) +
  geom_histogram(stat="bin", binwidth=1)+xlim(0,350)

max(mol2[,"num_reads"])
head(mol)
