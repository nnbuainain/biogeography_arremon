library("laser")### Get Branch Times from a newick tree
setwd("/Users/nnbuainain/Dropbox/Ornitologia_Nelson/Doutorado/UCE/Arremon/laser_package/may_2020/")
bt<-getBtimes("all.tre")
bt<-getBtimes("lowlandtree.tre")
bt<-getBtimes("uplandtree.tre")
bt<-getBtimes("tortree.tre")
bt<-getBtimes("bruntree.tre")
bt<-getBtimes("lystree.tre")

up<-read.tree("uplandtree.tre")
up<-force.ultrametric(up,method="nnls")
plot(up,cex=0.4)
write.tree(up,file = "uplandtree.tre")

## If somehow you ultrametric tree is shown as non ultrametric than you
##need to adjust the numeric precision

library(phytools)
tree<-read.tree("./arremon_uce_pl.tre")
ult<-force.ultrametric(tree) ##Default method nnls
write.tree(tree,"./ultrametric_to_laser_uce.tre",tree.names=TRUE)
bt<-getBtimes("./ultrametric_to_laser_uce.tre")

###Gamma statistics -If significant, we reject the null hypothesis 
##in favour of decreasing diversification rates through time.
###When gamma takes high positive values, one can infer an 
##increasing rate of diversification towards the present, 
##whereas high negative values indicate a slow- down of the diversification rate.

gamStat(bt)

### Lineages Through Time

plotLtt(bt)

### Test the fit of the data to various models with constant or variable rates through time

results<-fitdAICrc(bt, modelset = c("pureBirth", "bd","DDX", "DDL", "yule2rate"),ints="NULL")

###Fits SPVAR MODEL - fitSPVAR fits a model with an exponentially declining 
#speciation rate through time and constant extinction.

fitSPVAR(bt)

###Fits BOTH Models - fitBOTHVAR fits a model where both speciation and 
##extinction rates can vary through time.

fitBOTHVAR(bt)

#### Fits EXVAR Model - fitEXVAR fits a model with exponentially increasing
#extinction and constant speciation

fitEXVAR(bt)

## In case AIC is too close between two models you can calculate the weighted AIC
library(MuMIn)
#create a vector with aic values
aic<-c(-37.41,-35.430,-34.580)
#calculate aic
Weights(aic)
#divide greater value by second greater and have a proportion of 
#probabilities. if  aic1/aic2=2 the probility of model 1 being right is two
#times the probability of model 2.

