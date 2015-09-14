rm(list=ls())
source("SegFunction.R")
### load signal of WBS simulation
source("Simulation.R")
nbsimu <- length(Simu)
########################################################################
########################################################################
### 
set.seed(1)
e_coef <- 1 ### increase the size of the profile (increase the variance)

perfTrue <- list()


im <- 1
seg.method <- c("WBS.ssic", "WBS.bic", "WBS.mbic",
	"Fpop.3", "Fpop.2",  "Fpop.1",
	"Smuce.0.55", "Smuce.0.45", "Smuce.0.35")


getPerf <- function(i, nbrep=100){
  perfT <- list()
  for(approach in seg.method){
    perfT[[approach]] <- list()
  }

  Ktrue <- Simu[[i]]$Ktrue
  bkptrue <- as.integer(Simu[[i]]$bkpPage29[-c(1, Ktrue+1)])*e_coef
  signaltrue <-  rep(Simu[[i]]$signal, each=e_coef)
  sigmatrue <- Simu[[i]]$sigma *sqrt(e_coef)

  for(j in 1:nbrep){
    x <- signaltrue + rnorm(length(signaltrue), sd=sigmatrue)

    for(approach in seg.method){
     res <- assess_K_and_MSE(x, approach=approach, signaltrue, bkptrue, Ktrue, sigmatrue)
     perfT[[approach]][[j]] <- res
    }
  }

  for(approach in seg.method){
    perfT[[approach]] <-  do.call("rbind", perfT[[approach]])
  }
return(perfT)
}


####################################################################
####################################################################
e_coef <- 1
perfTrue <- mclapply(1:nbsimu, FUN=getPerf, nbrep=500, mc.cores=3)

save(perfTrue, file="Res_PerfWBSandSmuceSimu.Rdata")

perfSummary <- list()
for(i in 1:6){
perfSummary[[i]] <- do.call("rbind", lapply(perfTrue[[i]], colMeans))
}
SimuName <- c("WBS_1", "WBS_2", "WBS_2'(sigma as in Smuce)", "WBS_3", "WBS_4", "WBS_5")
for(i in 1:6){
write(SimuName[i], file="Res_SimuWBSandSmuce.csv", append=TRUE)
write.table(perfSummary[[i]], file="Res_SimuWBSandSmuce.csv", sep=";", append=TRUE, row.names=T, col.names=T)
}


