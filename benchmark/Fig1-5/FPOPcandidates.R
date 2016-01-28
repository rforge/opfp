setwd("~/Fig1-5")

source("FPOP-plot.R")  
set.seed(124)
y<-c(rnorm(60,2,1.5),rnorm(60,0,1))
pdf("FPOPcandidates.pdf",width=5,height=5)
OPRplot(y,log(100),79,min(y),max(y))
dev.off()


