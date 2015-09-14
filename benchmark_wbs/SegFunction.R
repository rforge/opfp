
require(fpop)
require(wbs)
require(stepR)
require(breakpointError)

################################################################
################################################################
varDiff <- function(x, method='MAD'){
  n = length(x)
  if(method == "MAD"){
	return(mad(diff(x)/sqrt(2)))	
  }
  if(method == "HALL"){
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(x)
    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]   
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
}

################################################################
################################################################
MSEf <- function(x, bkp, signal){
  bkp <- c(0, bkp)
  lgSeg <- diff(bkp)
  segNb <- rep(1:length(lgSeg), lgSeg)
  smt <- rep(by(x, segNb, FUN=mean)[1:length(lgSeg)], lgSeg)
  return(mean((signal-smt)^2))
}


##

fit <- function(x, bkp){
  bkp <- c(0, bkp)
  lgSeg <- diff(bkp)
  segNb <- rep(1:length(lgSeg), lgSeg)
  smt <- rep(by(x, segNb, FUN=mean)[1:length(lgSeg)], lgSeg)
  return(smt)
}

################################################################
################################################################
assess_K_and_MSE <- function(x, approach, signaltrue, bkptrue, Ktrue, sigmatrue, Kmax=50){
### WBS
if(grepl("^WBS", approach)){
   penalty <- paste(gsub("WBS.", "", approach), ".penalty", sep="")
   res <- wbs(x)
   bkp_predicted <- changepoints(res, Kmax=Kmax)$cpt.ic[[penalty]]
   if(is.na(bkp_predicted[1])){
	bkp_predicted <- length(x)
   } else {
        bkp_predicted <- sort(c(bkp_predicted, length(x)))
   }
}
### Fpop
if(grepl("^Fpop", approach)){
   penalty <- as.integer(gsub("Fpop.", "", approach))
   res <- Fpop(x/varDiff(x, method='MAD'), lambda=penalty*log(length(x)))
   bkp_predicted <- res$t.est 
}
#### Smuce
if(grepl("^Smuce", approach)){
  penalty <- as.numeric(gsub("Smuce.", "", approach))
  res <- smuceR(x, 1:length(x), family="gauss", alpha=penalty)
  bkp_predicted <- res$rightEnd

}

   ## nb breaks
   K <- length(bkp_predicted)
   CorrectSel <- (K == Ktrue)
   ## MSE
   MSE <- MSEf(x, bkp_predicted, signaltrue)/(sigmatrue^2)
   ## error
   bkp_error <- breakpointError(bkp_predicted[-K], bkptrue, length(x))
   error_b <- errorDetails(bkp_predicted[-c(1,K)], bkptrue, length(x))
   FP_ber = sum(error_b$false.positive)
   FN_ber = sum(error_b$false.negative)
   exactTP <- length(intersect(bkp_predicted[-K], bkptrue))
   exactFP <- K - exactTP
return(c(Khat=CorrectSel, MSE=MSE, exactTP=exactTP, exactFP=exactFP, BkpE=bkp_error, BerFP=FP_ber, BerFN=FN_ber))
}


