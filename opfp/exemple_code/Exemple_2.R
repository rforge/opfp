#### Guillem Rigaill Feb 2014
#### Example 2 : run op + fp, pelt and compare the runtime for different scenario

###### Exemple 2
###### Runtime comparaison on simulations
require(fpop)
require(changepoint)
require(wbs)
checkRuntime <- function(x, lambda){

  t1 <- system.time(res1 <- Fpop(x, lambda=lambda))[3]
  t2 <- system.time(res2 <- fpsn(x, 100))[3]
  t3 <- system.time(res3 <- cpt.mean(x, penalty="Manual", pen=lambda, method="PELT"))[3]
  ## res<- cpt.mean(x, Q=10, method="BinSeg", penalty="None")
  t4<- system.time( res4 <- multiBinSeg(x, 100))[3]
  t5 <- system.time( res5 <- wbs(x))[3]
  
  resT <- c(length(x), t1, t2, t3, t4, t5)

  return(resT)
}

###### with no break
res <- list()


ns <- rep(20*2^c(5:13), 8); i <- 1;
for(n in ns){
  if(i %% 2 == 1) print(paste("N =", n, "Turn : ", i, "/", length(ns)))
  x <- rnorm(n)
  res[[length(res)+1]] <- checkRuntime(x, lambda=log(n))
  i <- i+1
}


allTimes <- matrix(unlist(res), ncol=6, byrow=T)
mTimes <- aggregate(allTimes[, -1], by=list(allTimes[, 1]), FUN=mean)
png("img/Runtime_no_bkp-log.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], , type="b", lwd=5, lty=1:4, pch=c(20, 22, 24, 25, 23), col=1:5,  xlab="log-length", ylab="log-time", 
log="xy", cex.lab=1.7)
legend("topleft", legend=c("OP+FP", "SN+FP", "OP+IP", "BinSeg", "WBinSeg"), pch=c(20, 22, 24, 25, 23), col=1:5, cex=3, lwd=4)
dev.off()

#png("img/Runtime_no_bkp-log_withoutfpop.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -c(1, 2)], , type="b", lwd=5, lty=1:4, pch=c( 22, 24, 25), col=2:4,  xlab="log-length", ylab="log-time", 
#log="xy", cex.lab=1.7)
#legend("topleft", legend=c("SN+FP", "OP+IP", "BinSeg"), pch=c( 22, 24, 25), col=2:4, cex=3, lwd=4)
#dev.off()

#png("img/Runtime_no_bkp_raw.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -1], , type="b", lwd=5, lty=1:4, pch=c(20, 22, 24, 25), col=1:4,  xlab="log-length", ylab="log-time", cex.lab=1.7)
#legend("topleft", legend=c("OP+FP", "SN+FP", "OP+IP", "BinSeg"), pch=c(20, 22, 24, 25), col=1:4, cex=3, lwd=4)
#dev.off()
save(mTimes, file="mTimes_no_bkp.Rdata")

######################################
######################################
##### with increasing breaks
require(fpop)
require(changepoint)
checkRuntime <- function(x, lambda, K){

  t1 <- system.time(res1 <- Fpop(x, lambda=lambda))[3]
  ### 
  if(K <= 300){
  t2 <- system.time(res2 <- fpsn(x, K))[3] ### a priori to long
  } else {t2 <- NA}
  t3 <- system.time(res3 <- cpt.mean(x, penalty="Manual", pen=lambda, method="PELT"))[3]
  t4 <- system.time( res4 <- multiBinSeg(x, K))[3]
  t5 <- system.time( res5 <- wbs(x))[3] 
  resT <- c(length(x), t1, t2, t3, t4, t5)
  return(resT)
}

res <- list()
ns <- rep(20*2^c(6:13), 8); i <- 1;
for(n in ns){
  if(i %% 2 == 1) print(paste("N =", n, "Turn : ", i, "/", length(ns)))
  x <- rnorm(n) + 3*rep(rep(0:1, n/(2^6*10)), rep(2^5*10*c(1, 1), n/(2^6*10)))
  K = sum(diff(3*rep(rep(0:1, n/(2^6*10)), rep(2^5*10*c(1, 1), n/(2^6*10))))!=0)
  res[[length(res)+1]] <- checkRuntime(x, lambda=2*log(n),K)
  i <- i+1
}

allTimes <- matrix(unlist(res), ncol=6, byrow=T)
mTimes <- aggregate(allTimes[, -1], by=list(allTimes[, 1]), FUN=mean)
png("img/Runtime_with_bkp-log.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], , type="b", lwd=5, lty=1:4, pch=c(20, 22, 24, 25, 23), col=1:5,  xlab="log-length", ylab="log-time", 
log="xy", cex.lab=1.7)
legend("topleft", legend=c("OP+FP", "SN+FP", "OP+IP", "BinSeg", "WBinSeg"), pch=c(20, 22, 24, 25, 23), col=1:5, cex=3, lwd=4)
dev.off()

#png("img/Runtime_with_bkp-log_withoutfpop.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -c(1, 2)], type="b", lty=1:3,lwd=4, pch=c( 22, 24, 25), col=2:4, xlab="log-length", ylab="log-time",
# log="xy", cex.lab=1.7)
#legend("topleft", legend=c("SN+FP", "OP+IP", "BinSeg"), pch=c(22, 24, 25), col=2:4, cex=3, lwd=4)
#dev.off()

#png("img/Runtime_with_bkp_raw.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -1], , type="b", lty=1:4,lwd=4, pch=c(20, 22, 24, 25), col=1:4,  xlab="log-length", ylab="log-time",
# cex.lab=1.7)
#legend("topleft", legend=c("OP+FP", "pDPA", "PELT", "BinSeg"), pch=c(20, 22, 24, 25), col=1:4, cex=3, lwd=4)
#dev.off()
#save(mTimes, file="mTimes_nover100_bkp.Rdata")


###### with no break only 
#checkRuntime <- function(x, lambda){

#  t1 <- system.time(res1 <- opfp(x, lambda=lambda))[3]
  #t2 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="PELT"))[3]
#  t3 <- system.time( res3 <- multiSeg_Binary_Mean(x, 500 ))[3] ### fixed number of breaks
  #t3 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="BinSeg"))[3]
#  resT <- c(length(x), t1, t3)

#  return(resT)
#}
#res <- list()


#ns <- rep(20*2^c(5:21), 4); i <- 1;
#for(n in ns){
#  print(paste("N =", n, "Turn : ", i, "/", length(ns)))
#  x <- rnorm(n)
#  res[[length(res)+1]] <- checkRuntime(x, lambda=2*log(n))
#  i <- i+1
#}


#mTimes <- matrix(unlist(res), ncol=3, byrow=T)
#png("img/Runtime_no_bkp_alone-log.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="log-lenght", ylab="log-runtime", log="xy")
#legend("topleft", legend=c("op+fp", "BinSeg"), pch=c(20, 24), col=1:3)
#dev.off()
#png("img/Runtime_no_bkp_alone_raw.png", width=800, height=800)
#matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="lenght", ylab="runtime")
#legend("topleft", legend=c("op+fp", "BinSeg"), pch=c(20, 24), col=1:3)
#dev.off()
#save(mTimes, file="mTimes_no_bkp_alone.Rdata")

