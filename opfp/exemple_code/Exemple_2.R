#### Guillem Rigaill Feb 2014
#### Example 2 : run op + fp, pelt and compare the runtime for different scenario

###### Exemple 2
###### Runtime comparaison on simulations
require(opfp)
require(changepoint)
checkRuntime <- function(x, lambda){

  t1 <- system.time(res1 <- opfp(x, lambda=lambda))[3]
  t2 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="PELT"))[3]
  #t3 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="BinSeg"))[3]
  resT <- c(length(x), t1, t2)

  return(resT)
}

###### with no break
res <- list()


ns <- rep(20*2^c(3:12), 2); i <- 1;
for(n in ns){
  print(paste("N =", n, "Turn : ", i, "/", length(ns)))
  x <- rnorm(n)
  res[[length(res)+1]] <- checkRuntime(x, lambda=2*log(n))
  i <- i+1
}


mTimes <- matrix(unlist(res), ncol=3, byrow=T)
png("img/Runtime_no_bkp-log.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="log-lenght", ylab="log-runtime", log="xy")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
png("img/Runtime_no_bkp_raw.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="lenght", ylab="runtime")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
save(mTimes, file="mTimes_no_bkp.Rdata")

######################################
######################################
##### with increasing breaks
res <- list()
ns <- rep(20*2^c(3:18), 2); i <- 1;
for(n in ns){
  if(i %% 10 == 1) print(paste("N =", n, "Turn : ", i, "/", length(ns)))
  x <- rnorm(n) + 3*rep(rep(0:1, n/(2^4*10)), rep(2^3*10*c(1, 1), n/(2^4*10)))
  res[[length(res)+1]] <- checkRuntime(x, lambda=2*log(n))
  i <- i+1
}
mTimes <- matrix(unlist(res), ncol=3, byrow=T)
png("img/Runtime_nover100_bkp-log.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="log-lenght", ylab="log-runtime", log="xy")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
png("img/Runtime_nover100_bkp_raw.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="lenght", ylab="runtime")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
save(mTimes, file="mTimes_nover100_bkp.Rdata")


###### with no break only 
checkRuntime <- function(x, lambda){

  t1 <- system.time(res1 <- opfp(x, lambda=lambda))[3]
  #t2 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="PELT"))[3]
  #t3 <- system.time(res2 <- cpt.mean(x, penalty="Manual", pen=lambda, method="BinSeg"))[3]
  resT <- c(length(x), t1)

  return(resT)
}
res <- list()


ns <- rep(20*2^c(3:20), 4); i <- 1;
for(n in ns){
  print(paste("N =", n, "Turn : ", i, "/", length(ns)))
  x <- rnorm(n)
  res[[length(res)+1]] <- checkRuntime(x, lambda=2*log(n))
  i <- i+1
}


mTimes <- matrix(unlist(res), ncol=2, byrow=T)
png("img/Runtime_no_bkp_alone-log.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="log-lenght", ylab="log-runtime", log="xy")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
png("img/Runtime_no_bkp_alone_raw.png", width=800, height=800)
matplot(mTimes[, 1], mTimes[, -1], pch=c(20, 22), xlab="lenght", ylab="runtime")
legend("topleft", legend=c("op+fp", "pelt", "BinSeg"), pch=c(20, 22, 23), col=1:3)
dev.off()
save(mTimes, file="mTimes_no_bkp_alone.Rdata")

