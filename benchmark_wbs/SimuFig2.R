rm(list=ls())
source("SegFunction.R")
simu <- function(n, iK){
if(iK > 1){
	bkp <- sort(sample(1:(n-1), iK-1))
	lg <- diff(c(0, bkp, n)) ### rep(n/K, K)
	signal <- rep(rep(c(0, 1), length(lg))[1:length(lg)], lg)  [1:n]
} else {
	signal <- rep(0, n)		
}
		
x <- signal+ rnorm(n,sd=0.5)

return(list(x=x, signaltrue=signal, Ktrue=iK, bkptrue=bkp, sigmatrue=0.5))
}


############################################################

getSimu <- function(n, iK, nbrep=10){
	Simu <- list()
	for(i in iK){
	Simu[[i]] <- lapply(rep(n, nbrep), simu, iK =i)
	}
	return(Simu)
}
##################################
set.seed(1)
n <- 200000

#seg.method <- c("WBS.ssic", "Fpop.2", "Smuce.0.45")
seg.method <- c("WBS.ssic", "WBS.bic", "WBS.mbic",
	"Fpop.3", "Fpop.2",  "Fpop.1",
	"Smuce.0.55", "Smuce.0.45", "Smuce.0.35")
NbBkp <- c(10, 50, 100, 500, 1000)
#NbBkp <- c(50, 100, 500)
Simu <- getSimu(n, NbBkp, nbrep=50)

perfT <- list()
for(i in NbBkp) perfT[[i]] <- list()

for(i in NbBkp){
  print(paste("Nb Bkp", i))
  for(approach in seg.method){
   print(approach)
    system.time(perfTMP <- mclapply(Simu[[i]], FUN=function(x){
		return(assess_K_and_MSE(x$x, approach=approach, x$signaltrue, x$bkptrue, 
			x$Ktrue, x$sigmatrue, Kmax = max(i*6/5, 50) ))}, mc.cores=3))
   perfT[[i]][[approach]] <- do.call("rbind", perfTMP)
  }

}
save(perfT, file="Res_PerfSpeedSpeed.RData")


perfSummary <- list()
for(i in NbBkp){
perfSummary[[i]] <- do.call("rbind", lapply(perfT[[i]], colMeans))
}
SimuName <- paste("n= 200000, Bkp=", NbBkp)
ii <- 1
for(i in NbBkp){
write(SimuName[ii], file="Res_Simu_Speed.csv", append=TRUE)
ii <- ii+1
write.table(perfSummary[[i]], file="Res_Simu_Speed.csv", sep=";", append=TRUE, row.names=T, col.names=T)
}

