####Function for PELT - requires cost function####
PELTnoCP<-function(y,cost,pen,K=0){
  n<-length(y)
  lenR<-c()
  F<-c()
  F[1]<-(-pen)
  CP<-list()
  CP[[1]]<-0
  R<-list()
  R[[1]]<-0
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
  for (taustar in 1:n){
    temp<-c()
    for (tau in R[[taustar]]){
      temp[tau+1]<-F[tau+1]+cost(y,tau+1,taustar,S1,S2,SS1,SS2)+pen
    }
    F[taustar+1]<-min(temp,na.rm=T)
    taudash<-which(temp==min(temp,na.rm=T))-1
    CP[[taustar+1]]<-c(CP[[taudash+1]],taudash)
    temp2<-R[[taustar]]
    temp2<-temp2[(F[temp2+1]+cost(y,temp2+1,taustar,S1,S2,SS1,SS2)+K)<=F[taustar+1]]
    R[[taustar+1]]<-c(temp2,taustar)
    lenR[taustar]<-length(R[[taustar]])
  }
  return(lenR)
}


####Code for running the function####
#y<-SimChange(1000,c(1,3),1.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")


#PE<-PELTnoCP(y[[2]],negloglike.norm.mean2,log(1000),0);PE



####Cost functions####
negloglike.norm.mean2<-function(y,i,j,S1,S2,SS1,SS2,var=1.5){
  n<-length(y)
  len<-j-i+1
  S<-S1[n+1]-S2[n+1-j]-S1[i]
  SS<-SS1[n+1]-SS2[n+1-j]-SS1[i]
  return((len/2)*log(2*pi*var)+(1/(2*var))*(SS-(1/len)*S^2))
}

negloglike.norm.var2<-function(y,i,j,S1,S2,SS1,SS2,mean=0){
  n<-length(y)
  len<-j-i+1
  S<-S1[n+1]-S2[n+1-j]-S1[i]
  SS<-SS1[n+1]-SS2[n+1-j]-SS1[i]
  var<-0.0001
  if (len!=1){
    var<-(1/len)*SS-((1/len)*S)^2
  }
  return((len/2)*log(2*pi*var)+(1/(2*var))*(SS-2*mean*S+mean^2))
}