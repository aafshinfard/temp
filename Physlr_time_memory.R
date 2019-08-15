# Tank Gokce Dilek for the initial version of the script

library("ggplot2")
library("dplyr")
data<-read.csv("/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.time.min.tsv", sep = "\t" , header=T)

levels(data[,1]) = c(levels(data[,1]),"molecules_fast")
data[10,1]="molecules_fast"
#data<- data[seq(dim(df)[1],1),]
df <- data.frame(
  steps=data$Physlrstep,
  time=(data$Begin + data$End)/2,
  #group= c("cpp","cpp","cpp","cpp","cpp","cpp","py","py","py"),
  upper=data$End,
  lower=data$Begin
)

#pl = ggplot(df, aes(x=factor(steps), y=time))+ geom_linerange(aes(ymin=df$lower, ymax=df$upper))+ coord_flip() + theme_bw()+ geom_text(label=sprintf("%f min",data$Elapsedtime),, size=3, hjust=0, vjust=0, check_overlap = T )

level_order<-data$Physlrstep
#p <- ggplot(df, aes(x=factor(steps, level=level_order), y=time)) + geom_linerange(aes(ymin=lower, ymax=upper))+ coord_flip() + theme_bw() + ggtitle("Time Usage")+ geom_text(label=sprintf("%f min",data$Elapsedtime), size=3, hjust=0, vjust=0, check_overlap = T) + ylim(0,8)
p <- ggplot(df, aes(x=factor(steps), y=time))
myplot <- p+ geom_linerange(aes(ymin=df$lower, ymax=df$upper), size=4)+ coord_flip() + theme_bw() + ggtitle("Time Usage")+ geom_text(label=sprintf("%f min",data$Elapsedtime), size=4, hjust=0, vjust=-0.7, check_overlap = T) #+ ylim(0,8)
#png("/path/file", units="in",width=XX, height=XX, res=500)
myplot
#dev.off()
#p <- ggplot(df, aes(steps, time, colour= group))

#memory usage:
data2<-read.csv("/projects/btl_scratch/aafshinfard/projects/physlr2/physlr/data/hg004.mem.tsv", sep = "\t" , header=T)
levels(data2[,1]) = c(levels(data2[,1]),"molecules_fast")
data2[1,1]="molecules_fast"
df2 <- data.frame(
  steps=data2$Physlrstep,
  memory=(data2$Memory)/(1024*1024),
  #group=c("cpp","cpp","py","py","py"),
  lower=0,
  upper=(data2$Memory)/(1024*1024)
)
level_order<-data$Physlrstep
p <- ggplot(df2, aes(factor(steps, level=level_order),memory))
p+ geom_linerange(aes(ymin=lower, ymax=upper), size=3) +coord_flip()+ theme_bw()+ ggtitle("Memory Usage (MB)")+ geom_text(label=round((data2$Memory)/(1024*1024), digits=1), size=4, hjust=0.5, vjust=-0.8, check_overlap = T) + labs(y="Memory (GB)", x="Physlr step")


#fixed versions after investigating Lauren's scripts:

#time #2:
f1chr4<-read.csv("/home/gdilek/benchmarking/f1chr4.time.txt.minutes.tsv", sep = "\t" , header=T)
totaltime<-f1chr4 %>% summarize(totaltime = sum(Elapsedtime)) %>% 
  as.double() %>% round(digits=2)
f1chr4$Physlrstep<-factor(f1chr4$Physlrstep, levels=f1chr4$Physlrstep)
ggplot(f1chr4, aes(f1chr4$Elapsedtime, f1chr4$Physlrstep)) + geom_segment(size=2, aes(x = f1chr4$Begin, y = f1chr4$Physlrstep, xend = f1chr4$End, yend = f1chr4$Physlrstep))+ xlab("Time (min)") + ylab("Physlr step") + theme_bw() + theme(text=element_text(size=10))+ ggtitle(paste("Physlr: f1chr4 -k32 -w32 -n100 -N1000\t\tTotal time: ", totaltime, " min")) + xlim(-1,8)+ geom_text(size=2.5, aes(x=(f1chr4$Begin + f1chr4$End)/2, y=f1chr4$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(f1chr4$Elapsedtime, digits = 2), " min", sep=""))

#mem #2:
f1chr4<-read.csv("/home/gdilek/benchmarking/f1chr4.mem.txt.tsv", sep = "\t", header=T)
totalmem<-f1chr4 %>% summarize(totaltime = sum(Memory/1000)) %>% 
  as.double() %>% round(digits=2)
f1chr4$Physlrstep<-factor(f1chr4$Physlrstep, levels=f1chr4$Physlrstep)
ggplot(f1chr4, aes(f1chr4$Memory/1000, f1chr4$Physlrstep)) + geom_segment(size=2, aes(x=0, y=f1chr4$Physlrstep, xend = f1chr4$Memory/1000, yend = f1chr4$Physlrstep)) + xlab("Memory (MB)") +ylab("Physlr step") + theme_bw()+theme(text=element_text(size=10))+ggtitle(paste("Physlr: f1chr4 -k32 -w32 -n100 -N1000\t\tTotal memory: ", totaltime, " MB")) + geom_text(size=2.5, aes(x=f1chr4$Memory/1000, y=f1chr4$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(f1chr4$Memory/1000, digits = 2), " MB", sep=""))+xlim(-400,4000)

#time for filter-barcodes:
f1chr2R<-read.csv("/home/gdilek/physlr_task/withtime/f1chr2R.time.ms.tsv", sep = "\t", header=T)
totaltime<-f1chr2R %>% summarize(totaltime = sum(Elapsedtime/1000)) %>% 
  as.double() %>% round(digits=2)
f1chr2R$FilterbxStep<-factor(f1chr2R$FilterbxStep, levels=f1chr2R$FilterbxStep)
ggplot(f1chr2R, aes(f1chr2R$Elapsedtime/1000, f1chr2R$FilterbxStep)) + geom_segment(size=2, aes(x = f1chr2R$Begin/1000, y = f1chr2R$FilterbxStep, xend = f1chr2R$End/1000, yend = f1chr2R$FilterbxStep))+ xlab("Time (s)") + ylab("Filter-barcodes step") + theme_bw() + theme(text=element_text(size=10))+ ggtitle(paste("Filter-barcodes: f1chr2R -n100 -N900\t\tTotal time: ", totaltime, " s")) + geom_text(size=2.5, aes(x=((f1chr2R$Begin + f1chr2R$End)/2)/1000, y=f1chr2R$FilterbxStep), hjust=0.5, vjust=-1, label=paste(round(f1chr2R$Elapsedtime/1000, digits = 2), " s", sep=""))


#try time with colors:
#install.packages("RColorBrewer")
library("RColorBrewer")
group<-c("cpp","cpp","py","py","py", "py", "py", "py", "py", "py")
levels(group)
f1chr4<-read.csv("/projects/btl/gdilek/benchmarking/f1chr4/f1chr4.time.txt.minutes.tsv", sep = "\t" , header=T)
totaltime<-f1chr4 %>% summarize(totaltime = sum(Elapsedtime)) %>% 
  as.double() %>% round(digits=2)
f1chr4$Physlrstep<-factor(f1chr4$Physlrstep, levels=f1chr4$Physlrstep)
ggplot(f1chr4, aes(f1chr4$Elapsedtime, f1chr4$Physlrstep, colour=factor(group))) +scale_colour_brewer(palette = "Set2", name="a") + geom_segment(size=2, aes(x = f1chr4$Begin, y = f1chr4$Physlrstep, xend = f1chr4$End, yend = f1chr4$Physlrstep))+ xlab("Time (min)") + ylab("Physlr step") + theme_bw() + theme(text=element_text(size=10))+ ggtitle(paste("Physlr: f1chr4 -k32 -w32 -n100 -N1000\t\tTotal time: ", totaltime, " min")) + xlim(-1,8)+ geom_text(size=2.5, aes(x=(f1chr4$Begin + f1chr4$End)/2, y=f1chr4$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(f1chr4$Elapsedtime, digits = 2), " min", sep=""))

ggplot(f1chr4, aes(f1chr4$Elapsedtime, f1chr4$Physlrstep)) + geom_segment(size=2, aes(x = f1chr4$Begin, y = f1chr4$Physlrstep, xend = f1chr4$End, yend = f1chr4$Physlrstep))+ xlab("Time (min)") + ylab("Physlr step") + theme_bw() + theme(text=element_text(size=10))+ ggtitle(paste("Physlr: f1chr4 -k32 -w32 -n100 -N1000\t\tTotal time: ", totaltime, " min")) + xlim(-1,8)+ geom_text(size=2.5, aes(x=(f1chr4$Begin + f1chr4$End)/2, y=f1chr4$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(f1chr4$Elapsedtime, digits = 2), " min", sep=""))
