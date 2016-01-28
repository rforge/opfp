####Function using Rigaill's pruning methods on optimal partitioning, for normal data with change in mean and variance of 1.5####
OptimalPartitioningRigNoCP<-function(y,cost,pen){
  n<-length(y)
  lenLOC<-c()
  F<-c()
  F[1]<--pen
  CP<-list()
  CP[[1]]<-0
  quad<-function(x,A,B,C){
    return(A*x^2+B*x+C)
  }
  LOC<-c()
  LOC[1]<-0
  Set<-list()
   D<-c(min(y),max(y))
  Set[[1]]<-D
  a<-c(0);b<-c(0);c<-c(F[1]+pen)
  for (taustar in 1:n){
    temp<-c()
    for (tau in 1:length(LOC)){
        vs<-LOC[tau]
        a[vs+1]<-a[vs+1]+1/3;b[vs+1]<-b[vs+1]-(2/3)*y[taustar];c[vs+1]<-c[vs+1]+(1/2)*log(3*pi)+(1/3)*(y[taustar])^2
        temp[tau]<-NA
        if(a[vs+1]==0){temp[tau]<-c[vs+1]}
        else{
          for(i in 1:(length(Set[[vs+1]])/2)){
            if(Set[[vs+1]][2*i-1]<=-b[vs+1]/(2*a[vs+1])&-b[vs+1]/(2*a[vs+1])<=Set[[vs+1]][2*i]){
              temp[tau]<-quad(-b[vs+1]/(2*a[vs+1]),a[vs+1],b[vs+1],c[vs+1])}
          }
          if(is.na(temp[tau])){temp[tau]<-min(quad(Set[[vs+1]],a[vs+1],b[vs+1],c[vs+1]))}
        }
    }
    F[taustar+1]<-min(temp,na.rm=T)
    taudash<-LOC[(which(temp==F[taustar+1]))]
    CP[[taustar+1]]<-c(CP[[taudash+1]],taudash)
    a[taustar+1]<-0;b[taustar+1]<-0;c[taustar+1]<-F[taustar+1]+pen
    Set[[taustar+1]]<-D
    ####THE PRUNING BIT####
    for(tau in LOC){
      if((b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen))<0){I<-NA}
      else{
        I<-c((-b[tau+1]-sqrt(b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen)))/(2*a[tau+1]),(-b[tau+1]+sqrt(b[tau+1]^2-4*a[tau+1]*(c[tau+1]-F[taustar+1]-pen)))/(2*a[tau+1]))
      }
      Set[[tau+1]]<-in1(Set[[tau+1]],I) #intersect function for continuous sets
      if (is.na(Set[[tau+1]][1])){
        LOC<-setdiff(LOC,tau) #discrete set uses the setdiff function
      }
      Set[[taustar+1]]<-setdiff1(Set[[taustar+1]],I) #continuous set uses my setdiff1 function
    }
    if(!is.na(Set[[taustar+1]][1])){
      LOC<-c(LOC,taustar)}
    lenLOC[taustar]<-length(LOC)
    ####PRUNING BIT ENDS####
  }
  return(lenLOC)
}


####Code to run the function####
#y<-SimChange(1000,c(1,2,3,0,5,4,6,5,1),1.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")

#OPR<-OptimalPartitioningRigNoCP(y[[2]],negloglike.norm.mean2,log(1000));OPR
#abline(v=OPR,col="blue")
