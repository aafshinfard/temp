
A <- matrix(runif(2500), 50, 50)
B = A * t(A)
sum(A>0.5)
z = z2 = z3 = z4 = {}

###### a single-component subgraph
# try1
alpha = 0.6
#rowSums(x)
for (i in 1:1000) {
  n=50
  
  x <- matrix(sample(1:100, n*n, replace=T), n,n) 
  x4 = x3 = x2 = x
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind]
  x2[ind] <- t(x2)[ind] 
  x3[ind] <- t(x3)[ind] 
  
  x[x < (1-alpha)*100] = 0
  x[x >= (1-alpha)*100] = 1
  y=cosine(x)
  z=c(z,y)
  #hist(y,xlim=c(0,1),breaks = 20,col="blue")
  
  x2[x2 < (1-alpha)*100] = 0
  x2[x2 >= (1-alpha)*100] = 1
  diag(x2)=1
  y2=cosine(x2)
  z2=c(z2,y2)
  #hist(y2,xlim=c(0,1),breaks = 20,col="blue")
  
  x3[x3 < (1-alpha)*100] = 0
  y3=cosine(x3)
  z3=c(z3,y3)
  #hist(y3,xlim=c(0,1),breaks = 20,col="blue")
  
  x4[x4 < (1-alpha)*100] = 0
  diag(x4) = round(rowSums(x4)/(n-1))
  y4=cosine(x4)
  z4=c(z4,y4)
  #hist(y4,xlim=c(0,1),breaks = 20,col="blue")
}

hist(z,xlim=c(0,1),breaks = 40,col="blue")
hist(z2,xlim=c(0,1),breaks = 40,col="blue")
hist(z3,xlim=c(0,1),breaks = 40,col="blue")
hist(z4,xlim=c(0,1),breaks = 40,col="blue")



z={}
for (i in 1:1000) {
  n=50
  alpha = 0.4
  x <- matrix(sample(1:100, n*n, replace=T), n,n) 
  x4 = x3 = x2 = x
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind]
  
  x[x < (1-alpha)*100] = 0
  x[x >= (1-alpha)*100] = 1
  diag(x)=1
  #hist(y)
  
  x[1:round(n/2),(round(n/2)+1):(dim(x)[2])] = 0
  x[(round(n/2)+1):(dim(x)[2]),1:round(n/2)] = 0
  y=cosine(x)
  z=c(z,y)
  #hist(y,xlim=c(0,1),breaks = 20,col="blue")
  
}
z1 = z
hist(z1,xlim=c(0,1),breaks = 40,col="blue")
hist(z1[z1>0.0001],xlim=c(0,1),breaks = 40,col="blue")
z={}
noise = 2
for (i in 1:1000) {
  
  n=50
  #alpha = 0.4
  x <- matrix(sample(1:100, n*n, replace=T), n,n) 
  x4 = x3 = x2 = x
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind]
  
  x[x < (1-alpha)*100] = 0
  x[x >= (1-alpha)*100] = 1
  diag(x)=1
  #hist(y)
  xp = x
  x[1:round(n/2),(round(n/2)+1):(dim(x)[2])] = 0
  x[(round(n/2)+1):(dim(x)[2]),1:round(n/2)] = 0
  s1 = sample(1:round(n/2) , noise)
  s2 = sample((round(n/2)+1):(dim(x)[2]) , noise)
  x[s1,s2]=xp[s1,s2]
  x[s2,s1]=xp[s2,s1]
  
  y=cosine(x)
  z=c(z,y)
  #hist(y,xlim=c(0,1),breaks = 20,col="blue")
  
}
hist(z,xlim=c(0,1),breaks = 40,col="blue")
hist(z[z>0.001],xlim=c(0,1),breaks = 40,col="blue")
hist(y,xlim=c(0,1),breaks = 40,col="blue")
hist(y[y>0.001],xlim=c(0,1),breaks = 40,col="blue")

sum(y==0)
dim(y)


hist(z,xlim=c(0,1),breaks = 40,col="blue")
hist(z[z>0.001],xlim=c(0,1),breaks = 40,col="blue")
x

######################## interleaving node
z={}
z2={}
z3={}
z4={}
z5={}
z_1c={}
noiseCosine2={}
noiseCosine2_more={}
noiseCosine3={}
noiseCosine3_more={}
noiseCosine4_1={}
noiseCosine4_1_more={}
noiseCosine4_2={}
noiseCosine4_2_more={}
noiseCosine5_1={}
noiseCosine5_1_more={}
noiseCosine5_2={}
noiseCosine5_2_more={}
noise = 5
node_noise = 5
alpha = 0.4
n=50
trials = 10000
for (i in 1:trials) {
  ### building the graph itself
  print(i)
  x <- matrix(sample(1:100, n*n, replace=T), n,n) 
  x4 = x3 = x2 = x
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind]
  ### binary edges
  x[x < (1-alpha)*100] = 0
  x[x >= (1-alpha)*100] = 1
  ### set diag to 1 (self loop)
  diag(x)=1
  x_1c = x
  y_1c=cosine(x)
  z_1c = c(z_1c,y_1c)
  
  x_prime = x
  ### Building 2 components out of it (balanced)
  x[1:round(n/2),(round(n/2)+1):(dim(x)[2])] = 0
  x[(round(n/2)+1):(dim(x)[2]),1:round(n/2)] = 0
  x[round(n/2),] = 0
  x[,round(n/2)] = 0 
  y=cosine(x)
  z=c(z,y)
  
  x_prime = x
  ### add noise
  s1 = sample(1:(round(n/2)-1) ,noise ,replace = TRUE)
  s2 = sample((round(n/2)+1):(dim(x)[2]) ,noise ,replace = TRUE )
  x[(s1-1)*n + s2]= 1
  x[(s2-1)*n + s1]= 1
  x2 = x
  y2=cosine(x2)
  z2=c(z2,y2)
  offDiag2 = x2
  offDiag2[1:(round(n/2)-1),1:(round(n/2)-1)] = 0
  offDiag2[(round(n/2)+1):(dim(x2)[2]),(round(n/2)+1):(dim(x2)[2])] = 0
  noiseCosine2 = c( noiseCosine2 , y2[ offDiag2[,] > 0 ])
  #colSums(offDiag2[(round(n/2)+1):(dim(x2)[2]),1:(round(n/2)-1)]) Zarbedar y2[(round(n/2)+1):(dim(x2)[2]),1:(round(n/2)-1)] # no! need more checks
  for(i in 1:(round(n/2)-1) )
    if(sum(offDiag2[(round(n/2)+1):(dim(x2)[2]),i]) > 0)
      noiseCosine2_more = c(noiseCosine2_more, y2[(round(n/2)+1):(dim(x2)[2]),i] )
  for(i in (round(n/2)+1):(dim(x2)[2]) )
    if(sum(offDiag2[1:(round(n/2)-1),i]) > 0)
      noiseCosine2_more = c(noiseCosine2_more, y2[1:(round(n/2)-1),i] )
  #noiseCosine2_more = c(noiseCosine2_more, )
  
  temp = x_prime
  x_prime = x
  x = temp
  ### add interleaving nodes
  x[round(n/2),] = 0
  s1 = sample(1:(round(n/2)-1) , node_noise)
  x[round(n/2),s1] = 1
  x[s1,round(n/2)] = 1
  s2 = sample((round(n/2)+1):(dim(x)[2]) , node_noise)
  x[round(n/2),s2] = 1
  x[s2,round(n/2)] = 1
  x3 = x
  y3=cosine(x3)
  z3=c(z3,y3)
  noiseCosine3 = c(noiseCosine3,y3[25,x3[25,]>0])
  noiseCosine3_more = c(noiseCosine3_more,y3[25,])
  
  x = x_prime
  ### add interleaving nodes + noise
  x[round(n/2),] = 0
  s1 = sample(1:round(n/2) , node_noise)
  x[round(n/2),s1] = 1
  x[s1,round(n/2)] = 1
  s2 = sample((round(n/2)+1):(dim(x)[2]) , node_noise)
  x[round(n/2),s2] = 1
  x[s2,round(n/2)] = 1
  x4 = x
  y4=cosine(x4)
  noiseCosine4_1 = c(noiseCosine4_1,y4[25,x4[25,]>0])
  noiseCosine4_1_more = c(noiseCosine4_1_more,y4[25,])
  noiseCosine4_2 = c( noiseCosine4_2 , y4[ offDiag2[,] > 0 ])
  z4=c(z4,y4)
  for(i in 1:(round(n/2)-1) )
    if(sum(offDiag2[(round(n/2)+1):(dim(x)[2]),i]) > 0)
      noiseCosine4_2_more = c(noiseCosine4_2_more, y4[(round(n/2)+1):(dim(x2)[2]),i] )
  for(i in (round(n/2)+1):(dim(x)[2]) )
    if(sum(offDiag2[1:(round(n/2)-1),i]) > 0)
      noiseCosine4_2_more = c(noiseCosine4_2_more, y4[1:(round(n/2)-1),i] )
  x_prime = x
  ## connect neighbours of the neighbours:
  x2 = x
  for(i in 1:n){
    x2[i,colSums(x[x[i,]==1,])>1] = 1
  }
  x = x2
  x5 = x
  y5=cosine(x5)
  z5=c(z5,y5)
  noiseCosine5_1 = c(noiseCosine5_1,y5[25,x5[25,]>0])
  noiseCosine5_1_more = c(noiseCosine5_1_more,y5[25,])
  noiseCosine5_2 = c( noiseCosine5_2 , y5[ offDiag2[,] > 0 ])
  for(i in 1:(round(n/2)-1) )
    if(sum(offDiag2[(round(n/2)+1):(dim(x)[2]),i]) > 0)
      noiseCosine5_2_more = c(noiseCosine5_2_more, y5[(round(n/2)+1):(dim(x2)[2]),i] )
  for(i in (round(n/2)+1):(dim(x)[2]) )
    if(sum(offDiag2[1:(round(n/2)-1),i]) > 0)
      noiseCosine5_2_more = c(noiseCosine5_2_more, y5[1:(round(n/2)-1),i] )
  #hist(y,xlim=c(0,1),breaks = 20,col="blue")
}
hist(z_1c,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 1 component - ",noise,"",node_noise,"",alpha))
hist(z_1c[z_1c>0.001],xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 1 component - zoomed -",noise,"",node_noise,"",alpha))

#ggplot(data = as.data.frame(z_1c[]))+geom_histogram(binwidth = 0.01, mapping = aes((z_1c[])),alpha = 0.5)+geom_histogram(binwidth = 0.01, mapping = aes((y_1c[25,])),alpha = 0.5,color = "blue")

hist(z,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 components - ",noise,"",node_noise,"",alpha))
hist(z[z>0.001],xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 components - zoomed",noise,"",node_noise,"",alpha))
dat <- data.frame(cond = factor(c(rep("All", each=length(z) ))), value = c(z))
ggplot(dat[dat[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.8, binwidth=0.05, position="identity" )


hist(z2,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (edges)",noise,"",node_noise,"",alpha))
hist(z2[z2>0.001],xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (edges) - zoomed",noise,"",node_noise,"",alpha))
#dat2 <- data.frame(cond = factor(c(rep("All", each=length(z2) ),rep("noise", each=length(noiseCosine2)))), value = c(z2,noiseCosine2))
#ggplot(dat2[dat2[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.5, binwidth=0.01, position="identity" )
dat2 <- data.frame(cond = factor(c(rep("All", each=length(z2) ),rep("noise2", each=length(noiseCosine2_more)),rep("noise", each=length(noiseCosine2)))), value = c(z2,noiseCosine2_more,noiseCosine2))
ggplot(dat2[dat2[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.8, binwidth=0.05, position="identity" )

hist(z3,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes)",noise,"",node_noise,"",alpha))
hist(z3[z3>0.001],xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes) - zoomed",noise,"",node_noise,"",alpha))
#dat3 <- data.frame(cond = factor(c(rep("All", each=length(z3) ),rep("noise", each=length(noiseCosine3)))), value = c(z3,noiseCosine3))
#ggplot(dat3[dat3[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.5, binwidth=density(dat3$value)$bw, position="identity" )
dat3 <- data.frame(cond = factor(c(rep("All", each=length(z3) ),rep("noise2", each=length(noiseCosine3_more)),rep("noise", each=length(noiseCosine3)))), value = c(z3,noiseCosine3_more,noiseCosine3))
ggplot(dat3[dat3[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.8, binwidth=density(dat3$value)$bw * 2, position="identity" )

#ggplot(dat3, aes(x=value)) + geom_histogram(aes(y = ..density..), binwidth=density(dat3$value)$bw )
#ggplot(dat3[dat3[,"cond"]=="All" & dat3[,"value"]> 0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.5, binwidth=density(dat3$value)$bw )

#hist(z4,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes+edges)",noise,"",node_noise,"",alpha))
#hist(z4[z4>0.001],xlim=c(0,1),breaks = 40,col="blue",main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes+edges) - zoomed",noise,"",node_noise,"",alpha))
dat4 <- data.frame(cond = factor(c(rep("All", each=length(z4) ),rep("noise_more", each=length(noiseCosine4_1_more)),rep("noise", each=length(noiseCosine4_1)),rep("interl_Node_more", each=length(noiseCosine4_2_more)),rep("interl_Node", each=length(noiseCosine4_2)))), value = c(z4,noiseCosine4_1_more,noiseCosine4_1,noiseCosine4_2_more,noiseCosine4_2))
ggplot(dat4[dat4[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.6, binwidth= 0.05, position="identity" )

#hist(z5,xlim=c(0,1),breaks = 40,col="blue", main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes+edges) + NeiNei idea",noise,"",node_noise,"",alpha))
#hist(z5[z5>0.001],xlim=c(0,1),breaks = 40,col="blue",main = paste("Histogram of values in Sim mat - subgraph of 2 component + noise (nodes+edges) ) + NeiNei idea- zoomed",noise,"",node_noise,"",alpha))
dat5 <- data.frame(cond = factor(c(rep("All", each=length(z5) ),rep("noise_more", each=length(noiseCosine5_1_more)),rep("noise", each=length(noiseCosine5_1)),rep("interl_Node_more", each=length(noiseCosine5_2_more)),rep("interl_Node", each=length(noiseCosine5_2)))), value = c(z5,noiseCosine5_1_more,noiseCosine5_1,noiseCosine5_2_more,noiseCosine5_2))
ggplot(dat5[dat5[,"value"]>0 ,], aes(x=value, fill =cond)) + geom_histogram(alpha=0.6, binwidth= 0.05, position="identity" )



hist(z2,xlim=c(0,1),breaks = 40,col="blue")
hist(z2[z2>0.001],xlim=c(0,1),breaks = 40,col="blue")
hist(y,xlim=c(0,1),breaks = 40,col="blue")
hist(y[y>0.001],xlim=c(0,1),breaks = 40,col="blue")
