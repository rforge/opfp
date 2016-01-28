####### Code to generate Figure 4 from Maidstone et al. 2016, Comparison of the number of candidate #####
####### changepoints stored over time by FPOP and PELT. Averaged over 1000 data sets with changepoints ##
####### at t = 20; 40; 60 and 80. #######################################################################


setwd("~/Fig1-5")
source("FPOP-noCP.R") #version of fpop which stores number of candidate changes at each step
source("PELT-noCP.R") #version of pelt which stores number of candidate changes at each step

n<-1000
sumPELT<-0
sumOPR<-0
for(i in 1:n){
  y<-list()
  a<-c(-3,3,-3,3,-3)
  y[[2]]<-c(rnorm(20,a[1],1),rnorm(20,a[2],1),rnorm(20,a[3],1),rnorm(20,a[4],1),rnorm(20,a[5],1))
  y[[1]]<-seq(0,1000,by=200)
  
  PELT<-PELTnoCP(y[[2]],negloglike.norm.mean2,log(1000))
  sumPELT<-sumPELT+PELT
  
  OPR<-OptimalPartitioningRigNoCP(y[[2]],negloglike.norm.mean2,log(1000))
  sumOPR<-sumOPR+OPR
}
avePELT<-sumPELT/n
aveOPR<-sumOPR/n

#save(avePELT, aveOPR, file="ConsideredPELTandFPOP.RData")

library(ggplot2)
library(dplyr)
library(directlabels)

numcands<-data.frame(x = c(1:100, 1:100),y=c(aveOPR, avePELT), type = c(rep("FPOP",100),rep("PELT",100))) 

algo.colors <-
  c(
    pDPA="#1B9E77",
    PELT="#D95F02",
    FPOP="#7570B3",
    BinSeg="#E7298A",
    ##dnacopy.default="#66A61E",
    SNIP="#E6AB02",
    SMUCE="#A6761D", SMUCE.default="#A6761D",
    WBS="#666666", WBS.default="#666666")



p <- ggplot(numcands,aes(x = x, y = y, group = type)) + geom_line(aes(color = type))+
  theme_bw()+scale_color_manual(values=algo.colors)+scale_fill_manual(values=algo.colors)

#dl <- direct.label(p, "last.polygons")

#dl<-dl +xlab("Time") + ylab("Number of candidates being considered") +xlim(c(0,110))

p<-p +xlab("Time") + ylab("Number of candidates being considered") +xlim(c(0,110))

pdf("ConsideredPELTandFPOP.pdf",width=5,height=5)
print(p)
dev.off()
