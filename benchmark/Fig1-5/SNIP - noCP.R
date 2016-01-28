
###Segment Neighbourhood Search Code - requires cost function####
SNPnoCP<-function(y,cost=negloglike.norm.mean2,maxseg=4,kappa=0){
  n<-length(y)
  S1<-0
  for(i in 1:n){
    S1[i+1]<-S1[i]+y[i]
  }
  S2<-0
  for (i in 1:n){
    S2[i+1]<-S2[i]+y[n+1-i]
  }
  SS1<-0
  for(i in 1:n){
    SS1[i+1]<-SS1[i]+(y[i])^2
  }
  SS2<-0
  for (i in 1:n){
    SS2[i+1]<-SS2[i]+(y[n+1-i])^2
  }
  q1<-matrix(nrow=n,ncol=n)
  for (i in 1:(n)){
    for (j in (i):n){
      q1[i,j]<-cost(y,i,j,S1,S2,SS1,SS2)
    }
  }
  tau<-matrix(nrow=maxseg,ncol=maxseg-1)
  output<-matrix(nrow=maxseg,ncol=maxseg)
  output[1,1]<-q1[1,n]
  Q<-matrix(nrow=maxseg,ncol=n)
  Q[1,]<-q1[1,]
  metaR<-list()
  for (m in 2:maxseg){
    R<-list()
    R[[m]]<-m-1
    for (j in (m):n){
      temp<-c()
      for (v in R[[j]]){
        temp[v]<-Q[m-1,v]+q1[v+1,j]
      }
      Q[m,j]<-min(temp,na.rm=TRUE)
      temp2<-R[[j]]
      temp2<-temp2[Q[m-1,temp2]+q1[temp2+1,j]+kappa<=Q[m-1,j]]
      R[[j+1]]<-c(temp2,j)   
    }
    metaR[[m]]<-R
    tau[m,1]<-which(temp==Q[m,n])
    if (m>2){
      for (i in 2:(m-1)){
        temp<-c()
        for (v in metaR[[m-i+1]][[tau[m,i-1]]]){
          temp[v]<-Q[m-i,v]+q1[v+1,tau[m,i-1]]
        }
        tau[m,i]<-which(temp==min(temp,na.rm=T))
      }}
    output[m,1]<-Q[m,n]
  }
  output[1:maxseg,2:maxseg]<-tau
  lenR<-list()
  for (i in 2:maxseg){
    temp<-c()
    for (j in i:n){
      temp[j]<-length(metaR[[i]][[j]])      
    }
    lenR[[i]]<-temp
  }
  return(lenR)
}



####Cost function####
negloglike.norm.mean2<-function(y,i,j,S1,S2,SS1,SS2,var=1.5){
  n<-length(y)
  len<-j-i+1
  S<-S1[n+1]-S2[n+1-j]-S1[i]
  SS<-SS1[n+1]-SS2[n+1-j]-SS1[i]
  return((len/2)*log(2*pi*var)+(1/(2*var))*(SS-(1/len)*S^2))
}
