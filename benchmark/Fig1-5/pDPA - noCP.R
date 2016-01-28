####Function for running PDPA on Normal data with a change in mean and variance 1.5####
RigaillnoCP<-function(y,maxseg){
  n<-length(y)
  lenLOC<-list()
  S1<-0
  for(i in 1:n){
    S1[i+1]<-S1[i]+y[i]
  }
  SS1<-0
  for(i in 1:n){
    SS1[i+1]<-SS1[i]+(y[i])^2
  }
  C<-matrix(nrow=maxseg,ncol=n)
  for (t in 1:(n)){
    mu<-S1[t+1]/t
    C[1,t]<-t/2*log(3*pi)+1/3*(SS1[t+1]-2*mu*S1[t+1]+t*mu^2)
  }
  quad<-function(x,A,B,C){
    return(A*x^2+B*x+C)
  }
  tau<-matrix(nrow=maxseg,ncol=n)
  output<-matrix(nrow=maxseg,ncol=maxseg+1)
  output[1,1]<-C[1,n]
  LOC<-list()
  D<-c(min(y),max(y))
  for (m in 2:maxseg){
    Set<-list()
    LOC[[m]]<-m-1
    a<-c();b<-c();c<-c()
    a[m-1]<-0;b[m-1]<-0;c[m-1]<-C[m-1,m-1]  #min cost of m-1 CPs, last one at m-1
    Set[[m-1]]<-D
    templenLOC<-c()
    for (j in m:n){
      a[j]<-0;b[j]<-0;c[j]<-C[m-1,j]  #min cost of m-1 CPs, last one at j
      Set[[j]]<-D
      temp<-c()
      for (v in LOC[[m]]){
        a[v]<-a[v]+1/3;b[v]<-b[v]-(2/3)*y[j];c[v]<-c[v]+(1/2)*log(3*pi)+(1/3)*(y[j])^2
        if((b[v]^2-4*a[v]*(c[v]-C[m-1,j]))<0){I<-NA}
        else{
          I<-c((-b[v]-sqrt(b[v]^2-4*a[v]*(c[v]-C[m-1,j])))/(2*a[v]),(-b[v]+sqrt(b[v]^2-4*a[v]*(c[v]-C[m-1,j])))/(2*a[v]))
        }
        Set[[v]]<-in1(Set[[v]],I) #intersect function for continuous sets
        if (is.na(Set[[v]][1])){
          LOC[[m]]<-setdiff(LOC[[m]],v) #discrete set uses the setdiff function
        }
        Set[[j]]<-setdiff1(Set[[j]],I) #continuous set uses my setdiff1 function
      }
      if(!is.na(Set[[j]][1])){
        LOC[[m]]<-c(LOC[[m]],j)
      }
      temp<-c()
      for (v in 1:length(LOC[[m]])){
        vs<-LOC[[m]][v]
        temp[v]<-NA
        if(a[vs]==0){temp[v]<-c[vs]}
        else{
          for(i in 1:(length(Set[[vs]])/2)){
            if(Set[[vs]][2*i-1]<=-b[vs]/(2*a[vs])&-b[vs]/(2*a[vs])<=Set[[vs]][2*i]){
              temp[v]<--b[vs]/(2*a[vs])
              temp[v]<-quad(temp[v],a[vs],b[vs],c[vs])}
          }
          if(is.na(temp[v])){temp[v]<-min(quad(Set[[vs]],a[vs],b[vs],c[vs]))}
        }
      }
      C[m,j]<-min(temp,na.rm=T)
      tau[m,j]<-LOC[[m]][(which(temp==C[m,j]))]
      templenLOC[j]<-length(LOC[[m]])
    }
    lenLOC[[m]]<-templenLOC
    em<-m
    taustar<-c()
    taustar[m+1]<-n
    while(em >1){
      taustar[em]<-tau[em,taustar[em+1]]
      em<-em-1
    }
    output[m,1]<-C[m,n]
    output[m,2:(m+1)]<-taustar[-1]
  }
  return(lenLOC)
}


####Code to run Rigaill's method####
#y<-SimChange(1000,c(1,.4,.8,2,3,1.6,.3,.4,.2,1,4,.2,.8,1.2,1.4),.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")
#Q2<-Rigaill(y[[2]],14);Q2




####Functions for doing set operations on continuous sets, where a set is stored as a vector of endpoints####


####Union####
u1 <- function(A,B) 
{ 
  if(is.na(A[1])){
    return(B)
  }
  else if(is.na(B[1])){
    return(A)
  }
  else{
    Both<-c(A,B)
    m<-matrix(Both,ncol=2,byrow=T)
    m2<-m[order(m[,1]),]
    Both<-as.vector(t(m2)) #vector now of ordered pairs
    U<-c()
    U[1:2]<-Both[1:2]
    for(i in 2:((length(Both))/2)){
      if(U[length(U)-1]<=Both[2*i-1]&Both[2*i]<=U[length(U)]){
        #do nothing#
      }
      else if(U[length(U)-1]<=Both[2*i-1]&Both[2*i-1]<=U[length(U)]&U[length(U)]<=Both[2*i]){
        U[(length(U)-1):length(U)]<-c(U[length(U)-1],Both[2*i])
      }
      else{U[(length(U)+1):(length(U)+2)]<-Both[(2*i-1):(2*i)]}
    }
    return(U)}
} 

####Intersection####
in1<-function(A,B){
  if(is.na(A[1])|is.na(B[1])){return(NA)}
  else{
    U<-c()
    for(i in 1:(length(A)/2)){
      for(j in 1:(length(B)/2)){
        i1<-max(A[2*i-1],B[2*j-1])
        i2<-min(A[2*i],B[2*j])
        if(i1<=i2){U[(length(U)+1):(length(U)+2)]<-c(i1,i2)}
      }
    }
    if(is.null(U)){return(NA)}
    I<-u1(U,U[1:2])
    return(I)}
}

####Set Difference####
setdiff1<-function(A,B,M=10000){
  if(is.na(A[1])){return(NA)}
  else if(is.na(B[1])){return(A)}
  else{
    delta<-1/M
    m<-matrix(B,ncol=2,byrow=T)
    m2<-m[order(m[,1]),]
    B<-as.vector(t(m2)) #vector now of ordered pairs
    index<-rep(c(-delta,delta),length=length(B))
    BC<-c(-M,B+index,M)
    se<-in1(A,BC)
    return(se)}
}


