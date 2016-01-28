setwd("~/Fig1-5")
source("PDPA-plot.R")  
set.seed(124)
y<-c(rnorm(60,2,1.5),rnorm(60,0,1))
pdf("PDPAcandidates.pdf",width=5,height=5)
Rigaill(y,2,45,2,0.5,4.5)
dev.off()

