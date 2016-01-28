####### Code to generate Figure 5 from Maidstone et al. 2016, Comparison of the number of candidate #####
####### changepoints stored over time by pDPA and SNIP. Averaged over 1000 data sets with changepoints ##
####### at t = 20; 40; 60 and 80. #######################################################################


setwd("~/Fig1-5")
source("pDPA - noCP.R")
source("SNIP - noCP.R")

##4changesSegNeigh

n<-1000
sumrig<-list(0,0,0,0,0)
sumSNP<-list(0,0,0,0,0)
averigi<-list()
aveSNPi<-list()
averig<-list()
aveSNP<-list()
maxseg<-5
for(i in 1:n){
  y<-list()
  a<-c(-3,3,-3,3,-3)
  y[[2]]<-c(rnorm(20,a[1],1),rnorm(20,a[2],1),rnorm(20,a[3],1),rnorm(20,a[4],1),rnorm(20,a[5],1))
  y[[1]]<-seq(0,1000,by=200)
  
  riggy<-RigaillnoCP(y[[2]],maxseg)
  SNP<-SNPnoCP(y[[2]],maxseg=maxseg)
  for(j in 2:maxseg){
    sumrig[[j]]<-sumrig[[j]]+riggy[[j]]
    sumSNP[[j]]<-sumSNP[[j]]+SNP[[j]]
  }
}

for(j in 2:maxseg){
  averig[[j]]<-sumrig[[j]]/n
  aveSNP[[j]]<-sumSNP[[j]]/n
}

library(ggplot2)
library(dplyr)
library(directlabels)

maxseg=5
pdf("ConsideredSNIPandPDPA.pdf",width=3,height=4.5)
for(j in 2:maxseg){
  numcands<-data.frame(x = c(1:100, 1:100),y=c(aveSNP[[j]], averig[[j]]), type = c(rep("SNIP",100),rep("pDPA",100))) 
  
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
  
 # dl <- direct.label(p, "last.polygons")
  
  #dl<-dl +xlab("Time") + ylab("Number of candidates being considered") +xlim(c(0,130))+ylim(0,100)
  
  p<-p +xlab("Time") + ylab("Number of candidates being considered") +xlim(c(0,130))+ylim(0,100)
  
  print(p)
}
dev.off()



