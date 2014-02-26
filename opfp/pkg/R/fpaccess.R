## February 2014 Guillem Rigaill <rigaill@evry.inra.fr> 

retour_op <- function(path, i){
   chaine <- integer(1)
   chaine[1] <- length(path)
   j <- 2
   while(chaine[j-1] > 0){
	chaine[j] <- path[chaine[j-1]]
	j=j+1	
	}
   return(rev(chaine)[-1])
}

##### op function pruning access
##### return 
##### bkp position of the break
##### 
opfp <- function(x, lambda, mini=min(x), maxi=max(x)){
	n <- length(x)
    A <- .C("colibri_op_R_c", signal=as.double(x), n=as.integer(n), lambda=as.double(lambda),   min=as.double(mini), max=as.double(maxi), path=integer(n), cost=double(n)
	, PACKAGE="opfp")
    A$bkp <- retour_op(A$path, n)
    return(A);	
} 
##### nj function pruning access
##### return 
##### t.est a matrix of size Kmax x Kmax with the position of the best breaks for a segmentation in k from 1 to Kmax
##### J.est the cost of the best segmentation in k from 1 to Kmax
##### it is just a link to cghseg
njfp <- function(xe, Kmax){
	cghseg:::segmeanCO(x, Kmax)
}



