library(BAMMtools)
library("phytools")

setwd("./")
tree <- read.tree("all.tre")
edata <- getEventData(tree, eventdata = "event_data.txt", burnin=0.1)
tree = ladderize (tree, right = FALSE)

#Mean/Average model
plot.bammdata(edata, lwd=2.5, method="phylogram", pal="RdBu",spex="netdiv",legen=T)
axisPhylo()

#Posterior of each number of shift
shift_probs <- summary(edata)

#Compute BayesFactor to see what model (number of shifts) fit better the data
postfile<-"./mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.3)
bfmat
#"The 95% credible set is the set of distinct shift configurations that account
#for 95% of the probability of the data."

css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5,set.limit=0.95)

#how many distinct shifts
css$number.distinct

#summary
summary(css)

plot.credibleshiftset(css,lwd=2.5,pal="RdBu")
axisPhylo()
###General reccomendations are to present the bestshiftconfiguration
#This shift configuration is the one with the maximum a posteriori (MAP) probability.

best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot.bammdata(best, lwd=2.5,pal="RdBu",method = "phylogram",spex = "netdiv")
addBAMMshifts(best, cex=2)
axisPhylo()

###Get Cohort Matrix, not relevant to me
#cmat <- getCohortMatrix(edata)
#cohorts(cmat, edata, lwd=3, pal="temperature", use.plot.bammdata=TRUE)


#Get Rates through time
##First for all tree
#Open a window with three slots
plot.new()
par(mfrow=c(1,3))
#especiation
st <- max(branching.times(tree))
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73",start.time=st,ratetype = "speciation",ylim=c(0,1))
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4))
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",ylim=c(0,1))



###Plot and Anotate node to perform analysis for each clade of interest
plot.phylo(tree)
nodelabels(cex = .5, bg = "yellow")

##Lowland (Arremon)
#Open a window with three slots

plot.new()
par(mfrow=c(1,3))
#mtDNA=62
#uce=64
#especiation
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73", start.time=st,ratetype = "speciation",ylim=c(0,1),node=59)
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4),node=59)
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",node=59,ylim=c(0,1))


##upland (brunneinucha+torquatus+Lysurus)
#mtDNA=40
#uce=41
#Open a window with three slots
plot.new()
par(mfrow=c(1,3))
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73", start.time=st,ratetype = "speciation",ylim=c(0,1),node=39)
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4),node=39)
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",node=39,ylim=c(0,1))

## Only brunneinucha
#Open a window with three slots
#mtDNA=41
#uce=53
plot.new()
par(mfrow=c(1,3))
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73", start.time=st,ratetype = "speciation",ylim=c(0,1),node=52)
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4),node=52)
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",ylim=c(0,1),node=52)

## Only torquatus
#Open a window with three slots
#mtDNA=52
#uce=43
plot.new()
par(mfrow=c(1,3))
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73", start.time=st,ratetype = "speciation",ylim=c(0,1),node=43)
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4),node=43)
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",ylim=c(0,1),node=43)

## Only Lysurus
#Open a window with three slots
#mtdna=61
#uce=52
plot.new()
par(mfrow=c(1,3))
plotRateThroughTime(edata, intervalCol="#009E73", avgCol="#009E73", start.time=st,ratetype = "speciation",ylim=c(0,1),node=41)
#Extinction
plotRateThroughTime(edata, intervalCol="#0072B2", avgCol="#0072B2", start.time=st,ratetype = "extinction",ylim=c(0,.4),node=41)
#Net diversification
plotRateThroughTime(edata, intervalCol="black", avgCol="black", start.time=st,ratetype = "netdiv",ylim=c(0,1),node=41)



# Results for rates in general

allrates <- getCladeRates(edata)
mean(allrates$lambda)

quantile(allrates$lambda, c(0.05, 0.95))


#Plotar node labels
plot.phylo(tree)
nodelabels(cex = .75, bg = "yellow")

#Clade Specific rates
#lowland
low <- getCladeRates(edata, node=59)
mean(low$lambda)
quantile(low$lambda, c(0.05, 0.95))

#upland
up <- getCladeRates(edata, node=39)
mean(up$lambda)
quantile(up$lambda, c(0.05, 0.95))

#brunneinucha
brun <- getCladeRates(edata, node=52)
mean(brun$lambda)
quantile(brun$lambda, c(0.05, 0.95))

#torquatus
tor <- getCladeRates(edata, node=43)
mean(tor$lambda)
quantile(tor$lambda, c(0.05, 0.95))

#Lysurus
lys <- getCladeRates(edata, node=41)
mean(lys$lambda)
quantile(lys$lambda, c(0.05, 0.95))


#PlotTree with branch_priors
branch_priors <- getBranchShiftPriors(edata, expectedNumberOfShifts = 1)
plot(branch_priors)

data(edata, eventdata)
edata <- getEventData(edata, eventdata, burnin=0.1)
branch_priors <- getBranchShiftPriors(edata, expectedNumberOfShifts = 1)
mo <- marginalOddsRatioBranches(edata, branch_priors)


####Plot rates densityplot 
#creta one data frame for each categry
low<-data.frame(clade=c("low"),rate=low)
up<-data.frame(clade=c("up"),rate=up)
brun<-data.frame(clade=c("brun"),rate=brun)
tor<-data.frame(clade=c("tor"),rate=tor)
lys<-data.frame(clade=c("lys"),rate=lys)

#combine them
library(ggplot2)
lowuprates<-rbind(low,up)
allrates<-rbind(low,brun,tor,lys)

#diversification rate
ggplot(allrates) +
  aes(x = rate.lambda, fill = clade) +
  geom_density(adjust = 1L,alpha=0.3) +
  scale_fill_hue() +
  theme_minimal()+ scale_fill_manual(values=c("#009E73", "#E69F00", "#56B4E9", "#CC79A7"))+
  theme_classic()

ggplot(lowuprates) +
  aes(x = rate.lambda, fill = clade) +
  geom_density(adjust = 1L,alpha=0.3) +
  scale_fill_hue() +
  theme_minimal()+ scale_fill_grey() + theme_classic()
  
#extinction rate
ggplot(allrates) +
  aes(x = rate.mu, fill = clade) +
  geom_density(adjust = 1L,alpha=0.3) +
  scale_fill_hue() +
  theme_minimal()+ scale_fill_manual(values=c("#009E73", "#E69F00", "#56B4E9", "#CC79A7"))+
  theme_classic()

ggplot(lowuprates) +
  aes(x = rate.mu, fill = clade) +
  geom_density(adjust = 1L,alpha=0.3) +
  scale_fill_hue() +
  theme_minimal()+ scale_fill_grey() + theme_classic()
