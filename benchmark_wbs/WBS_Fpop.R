rm(list=ls())
source("SegFunction.R")
### load signal of WBS simulation
source("Simulation.R")
nbsimu <- length(Simu)
require(exactRankTests)
e_coef <- 1
SimuName <- c("WBS_1", "WBS_2 (sigma as in WBS)", "WBS_2'(sigma as in Smuce)", "WBS_3", "WBS_4", "WBS_5")
set.seed(1)
pdf("MSE_WBSssic_Fpop2 .pdf")

for(i in 1:6){
print(i)

fpopRes <- list()
wbsRes <- list()
nbrep <- 2000
signals <- list()
for(jj in 1:nbrep){
  Ktrue <- Simu[[i]]$Ktrue
  bkptrue <- as.integer(Simu[[i]]$bkpPage29[-c(1, Ktrue+1)])*e_coef
  signaltrue <-  rep(Simu[[i]]$signal, each=e_coef)
  sigmatrue <- Simu[[i]]$sigma *sqrt(e_coef)
  signals[[jj]] <- signaltrue + rnorm(length(signaltrue), sd=sigmatrue)
}

  fpopRes <- mclapply(signals, FUN=assess_K_and_MSE, approach="Fpop.2", signaltrue, bkptrue, Ktrue, sigmatrue, mc.cores=3)
  wbsRes <- mclapply(signals, FUN=assess_K_and_MSE, approach="WBS.ssic", signaltrue, bkptrue, Ktrue, sigmatrue, mc.cores=3)
  smuceRes <- mclapply(signals, FUN=assess_K_and_MSE, approach="Smuce.0.45", signaltrue, bkptrue, Ktrue, sigmatrue, mc.cores=3)


resFpop <- do.call("rbind", fpopRes)
resWBS <- do.call("rbind", wbsRes)
resSmuce <- do.call("rbind", smuceRes)

### WBS
pvalW <- signif(wilcox.exact(resWBS[, 2]-resFpop[, 2], alternative="greater")$p.value,3)
pvalT <- signif(t.test(resWBS[, 2]-resFpop[, 2], alternative="greater")$p.value, 3)
med <- signif(mean(resWBS[, 2]-resFpop[, 2]),3)
WBSwins <- sum(resWBS[, 2] < resFpop[, 2])/nbrep
WBSwinsOrEq <- sum(resWBS[, 2] <= resFpop[, 2])/nbrep
slope <- signif(coefficients(model <- lm(resFpop[, 2] ~ resWBS[, 2] -1)),3)

### PLOT ##########
plot(resWBS[, 2], resFpop[, 2], xlab="MSE.wbs.ssic", ylab="MSE.fpop.2", main=paste(SimuName[i], 
	", MSE"), sub=paste("Mean(Wbs-Fp):", med, ", Ww:", WBSwins, 
	" (<) /", WBSwinsOrEq, "(<=)", ", pv W/T:", pvalW, "/", pvalT), pch=20, col="blue"); 
abline(a=0, b=1, col="red", lwd=1, lty=2)
abline(a=0, b=slope, lwd=1, lty=1, col="red")
legend("topleft", c("y=x", paste("y=", slope, "x", sep="")), lty=c(2, 1), col="red")

###########################
### Smuce
pvalW <- signif(wilcox.exact(resSmuce[, 2]-resFpop[, 2], alternative="greater")$p.value,3)
pvalT <- signif(t.test(resSmuce[, 2]-resFpop[, 2], alternative="greater")$p.value, 3)
med <- signif(mean(resSmuce[, 2]-resFpop[, 2]),3)
WBSwins <- sum(resSmuce[, 2] < resFpop[, 2])/nbrep
WBSwinsOrEq <- sum(resSmuce[, 2] <= resFpop[, 2])/nbrep
slope <- signif(coefficients(model <- lm(resFpop[, 2] ~ resSmuce[, 2] -1)),3)

### PLOT ##########
plot(resSmuce[, 2], resFpop[, 2], xlab="MSE.smuce.0.45", ylab="MSE.fpop.2", main=paste(SimuName[i], 
	", MSE"), sub=paste("Mean(Smuce-Fp):", med, ", Sw:", WBSwins, 
	" (<) /", WBSwinsOrEq, "(<=)", ", pv W/T:", pvalW, "/", pvalT), pch=20, col="blue"); 
abline(a=0, b=1, col="red", lwd=1, lty=2)
abline(a=0, b=slope, lwd=1, lty=1, col="red")
legend("topleft", c("y=x", paste("y=", slope, "x", sep="")), lty=c(2, 1), col="red")
### END PLOT ######
}
dev.off()
