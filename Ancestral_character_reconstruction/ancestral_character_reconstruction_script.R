library(phytools)
library(ape)
library("geiger")
setwd("./")
tree<-read.tree("./data/all.tre")

x<-read.table("./data/habitat_3.txt",head=FALSE)
#or
x<-factor(c("montane","montane","montane","montane","montane","montane","montane","montane","montane","montane","montane_humidlowland","montane_humidlowland","montane","montane","montane","montane","montane","montane","montane","montane","montane","montane","montane","montane","humidlow","humidlow","humidlow","drylow","drylow","drylow","drylow","dry_humid_low","dry_humid_low","montane_humidlowland","drylow","drylow","drylow","drylow","drylow"))

ace<-ace(x,tree,type="discrete")
cols<-setNames(c("red","blue","green","yellow","black"),levels(x))
plot(tree,cex=.4)
tiplabels(pch = 21, bg = cols[as.numeric(x)], cex = 1, adj = 1)
nodelabels(pie = ace$lik.anc, piecol = cols, cex = 0.6)

nodelabels(node=1:tree$Nnode+Ntip(tree),pie=fitER$lik.anc,piecol=cols,cex=0.4)


###We can also consider a model in which the backward & forward rates between states
##are permitted to have different values:


fitARD<-ace(svl,tree,model="ARD",type="discrete")
fitARD
fitARD$lik.anc
plot(tree,cex=.4)
tiplabels(pch = 21, bg = cols[as.numeric(x)], cex = 1, adj = 1)
nodelabels(pie = fitARD$lik.anc, piecol = cols, cex = 0.6)

####ALTERNATIVE

svl<-read.csv("./habitat_data.csv",row.names=1,sep=";")
svl<-as.matrix(svl)[,2]
svl
mtree<-make.simmap(tree,svl,model="ER")
plot(mtree)

###### Testing Phylogenetic signal with lambda package geiger

fitDiscrete(tree, svl, model="ER", transform = "lambda", data.names=NULL, plotlnl=F, qLimits=c(0.0001, 1000), pLimits=c(0.00001, 10))

###### Testing Phylogenetic signal with D statistics package caper, needs data to be in binary format
## see Zhuo Chen1,2 & John J. Wiens 
## https://doi.org/10.1038/s41467-020-14356-3

?phylo.d()




#SELECTING THE BEST MODEL with the weighted AIC

ER<-fitMk(tree,svl,model = c("ER"))
ARD<-fitMk(tree,svl,model = c("ARD"))
SYM<-fitMk(tree,svl,model = c("SYM"))
aic<-setNames(c(AIC(ER),AIC(ARD),
                AIC(SYM)),c("ER","ARD","SYM"))
##Choose the model with the lowest AIC value
aic

#or the model with the AICw closest to 1
aic.w(aic)

#Adapt the object cols according to the dataset color attributed to levels 

y<-factor(c("dry","humid","humid_dry"))
cols<-setNames(c("#999999","#009E73","#E69F00"),levels(y))

y<-factor(c("ausente","intermediario","polimorfico","presente"))
cols<-setNames(c("white","#999999","#0072B2","black"),levels(y))

y<-factor(c("black","orange","orange_black","white_gray","yellow_black"))
cols=setNames(c("black","#E69F00","#D55E00","#999999","#F0E442"),levels(y))
y<-factor(c("green","gray"))
cols=setNames(c("#999999","#009E73"),levels(y))
y<-factor(c("absent","complete","reduced","reduced_absent","reduced_complete"))
cols=setNames(c("white", "black", "#999999", "#009E73", "#0072B2"),levels(y))
y<-factor(c("black","black or gray","white","white and gray","white and green"))
cols=setNames(c("black", "#0072B2", "white", "#999999", "#009E73"),levels(y))
y<-factor(c("absent","present"))
cols=setNames(c("white", "black"),levels(y))
y<-factor(c("no","yes","polimorphic"))
cols=setNames(c("white","#999999", "black"),levels(y))
y<-factor(c("black","gray","white"))
cols=setNames(c("black","#999999", "white"),levels(y))
y<-factor(c("no","yes"))
cols=setNames(c("white","black"),levels(y))
y<-factor(c("montane","lowland","montane_lowland"))
cols=setNames(c("#009E73","#0072B2","#CC79A7"),levels(y))


###

mtrees<-make.simmap(tree,svl,model="SYM",nsim=10000)
par(mfrow=c(10,10))
null<-sapply(mtrees,plotSimmap,colors=cols,lwd=1,ftype="off")
pd<-summary(mtrees)
pd
dev.off()

plot(pd,fsize=0.6,ftype="i",ylim=c(-2,Ntip(tree)))
add.simmap.legend(prompt=FALSE,x=0,y=0,vertical=FALSE,shape="circle",colors)


#write results to file

sink("./results/habitat_montane_lowland.txt",append = TRUE)
mtrees
pd
sink()

pdf("./results/habitat_humid_dry.pdf",width = 6,height = 4.63,useDingbats=FALSE)
plot(pd,fsize=0.6,ftype="i",ylim=c(-2,Ntip(tree)))
axis(1,pos=-0.05*h,lwd=2)
dev.off()

pdf("./results/habitat_montane_lowland_s.pdf",width = 4,height = 3,useDingbats=FALSE)
plot(pd,fsize=0.6,ftype="i",ylim=c(-2,Ntip(tree)))
axis(1,pos=-0.05*h,lwd=2)
dev.off()

pdf("./results/habitat_humid_dry_round.pdf",width = 6,height = 4.63,useDingbats=FALSE)
plot(pd,type="fan",fsize=0.6,ftype="i",lwd=1.5,cex=c(0.5,0.3))
dev.off()

#in Case you want to choose the colors
plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(tree)))
add.simmap.legend(colors=cols[2:1],prompt=FALSE,x=0,y=-4,vertical=FALSE)

plotTree(pd,type="fan",fsize=0.6,ftype="i",lwd=1)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitARD$lik.anc,piecol=cols,cex=0.4)

##Colors chosen and time scale

pdf("./results/humid_dry.pdf",width = 6,height = 4.63,useDingbats=FALSE)
plot(pd,fsize=0.6,ftype="i",colors=cols,ylim=c(-2,Ntip(tree)))
axis(1,pos=-0.05,lwd=1)
dev.off()

#Round
pdf("./results/humid_dry_round.pdf",width = 6,height = 4.63,useDingbats=FALSE)
plot(pd,type="fan",fsize=0.6,colors=cols,ftype="i",lwd=1.5,cex=c(0.5,0.3))
dev.off()

#plot with posterior along branches
pdf("./montane_humid_posterior.pdf",width = 6,height = 4.63,useDingbats=FALSE)
plotSimmap(mtrees[[1]],cols,fsize=0.6,ftype="i",ylim=c(-2,Ntip(tree)))
nodelabels(pie=pd$ace,piecol=cols,cex=0.5)
axis(1,pos=-0.05,lwd=1)
dev.off()
