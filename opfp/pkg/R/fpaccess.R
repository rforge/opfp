retour_op <- function
### This function is use by the fpop function to recover the best segmentation from 1:n from the C output
(path, 
### the path vector of the "colibri_op_R_c C" function
i
### the last position to consider in the path vector
){
   chaine <- integer(1)
   chaine[1] <- length(path)
   j <- 2
   while(chaine[j-1] > 0){
	chaine[j] <- path[chaine[j-1]]
	j=j+1	
	}
   return(rev(chaine)[-1])
### return a vector with the best change-points w.r.t. to L2 to go from point 1 to i
}

Fpop <- function
### Function calling the fpop algorithm, use functional pruning and optimal partionning. It is implemented for the L2-loss functon
(x, 
### A vector of double : the signal to be segmented
lambda, 
### Value of the penalty
mini=min(x), 
### Min value for the mean parameter of the segment
maxi=max(x)
### Max value for the mean parameter of the segment
){
  n <- length(x)
  A <- .C("colibri_op_R_c", signal=as.double(x), n=as.integer(n), 
		lambda=as.double(lambda),   min=as.double(mini), 
		max=as.double(maxi), path=integer(n), cost=double(n)
	, PACKAGE="fpop")
    A$bkp <- retour_op(A$path, n)
    return(A);	
### return a list with a vector containing the position of the change-points
} 

fpop_analysis <- function
### A function to count the number of intervals and or candidate segmentation at each step of fpop (under-developpemment)
(x,
### A vector of double : the signal to be segmented
lambda, 
### Value of the penalty
mini=min(x), 
### Min value for the mean parameter of the segment
maxi=max(x)
### Max value for the mean parameter of the segment
){

	n <- length(x)
    A <- .C("colibri_op_R_c_analysis", signal=as.double(x), n=as.integer(n), lambda=as.double(lambda),   min=as.double(mini), max=as.double(maxi), path=integer(n), cost=double(n), nbCandidate=integer(n)
	, PACKAGE="fpop")
    A$t.est <- retour_op(A$path, n)
    return(A);	
### return a list with a vector containing the position of the change-points t.est
} 

##### nj function pruning access
##### return 
##### t.est a matrix of size Kmax x Kmax with the position of the best breaks for a segmentation in k from 1 to Kmax
##### J.est the cost of the best segmentation in k from 1 to Kmax
##### it is just a link to cghseg
fpsn <- function
### Function to run the pDPA algorithm with the L2 loss (it is a wrapper to cghseg for now)
(x, 
### A vector of double : the signal to be segmented
Kmax
### Maximum number of segments. The snfp will recover all the best segmentation w.r.t. to the L2 in 1 to Kmax segments
){
	cghseg:::segmeanCO(x, Kmax)
### return a list with a J.est vector containing the L2 loss and a t.est matrix with the changes of the segmentations in 1 to Kmax
}



