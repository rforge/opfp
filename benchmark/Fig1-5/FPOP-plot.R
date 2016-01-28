
####Function using Rigaill's pruning methods on optimal partitioning, for normal data with change in mean and variance of 1.5####
OPRplot<-function(y,pen,timestep,mu.min,mu.max){
  n<-length(y)
  F<-c()
  F[1]<--pen
  CP<-list()
  CP[[1]]<-
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
    OLOC=c(LOC,taustar)
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
    ####PRUNING BIT ENDS####
    ######PLOT timestep-1 ###########
    if(taustar==timestep-1){
    x=seq(mu.min,mu.max,length=100)
    co=1:length(OLOC)
    OLOCdash<-c()
    for(it in OLOC){
      if(is.na(Set[[it+1]][1])==FALSE){
        OLOCdash<-c(OLOCdash,it)}
    }
    for(kk in OLOCdash){
      if(kk==OLOCdash[1]){
        plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[(kk)==OLOCdash],ylim=F[taustar+1]+c(0,10),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
      }else{
        lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[(kk)==OLOCdash])
      }
      if(length(Set[[kk+1]])>0){
        mm=length(Set[[kk+1]])/2
        Set2<-Set
        if(is.na(Set[[kk+1]][1])==FALSE){
          if(Set[[kk+1]][1]<mu.min){Set2[[kk+1]][1]=mu.min}
          if(Set[[kk+1]][length(Set[[kk+1]])]>mu.max){Set2[[kk+1]][length(Set[[kk+1]])]=mu.max}}
        delta=.4 #size of "fleck"
        delta2=.005 #size of space between "flecks"
        for(ii in 1:mm) {lines(Set2[[kk+1]][2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1])),2),lwd=3,col=co[(kk)==OLOCdash])
                         #browser()
                         lines(rep(Set2[[kk+1]][2*ii-1],2)+delta2,c(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))-delta,min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))+delta),lwd=3,col=co[(kk)==OLOCdash])
                         lines(rep(Set2[[kk+1]][2*ii-0],2)-delta2,c(min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))-delta,min(c[OLOCdash+1]-0.25*b[OLOCdash+1]^2/(1e-3+a[OLOCdash+1]))+delta),lwd=3,col=co[(kk)==OLOCdash])
                         ######boldbit
                         if(is.na(Set2[[kk+1]][2*ii-1])==FALSE){
                           x2<-seq(Set2[[kk+1]][2*ii-1],Set2[[kk+1]][2*ii-0],length=50)
                           lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[(kk)==OLOCdash],lwd=4)}
                         ######
        }
      }
    }
    legend("topleft", legend=paste(OLOCdash," "), col=co,pch=15,bg="white")}
  ####PLOT END########
  ######PLOT timestep mid###########
  if(taustar==timestep){
  x=seq(mu.min,mu.max,length=100)
  OLOCdash<-c(OLOCdash,taustar)
  co=1:length(OLOCdash)
  for(kk in OLOC){
    if(kk==OLOC[1]){
      plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==OLOC],ylim=F[taustar+1]+c(0,10),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
    }else{
      lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==OLOC])
    }
    if(length(Set[[kk+1]])>0){
      mm=length(Set[[kk+1]])/2
      Set2<-Set
      if(is.na(Set[[kk+1]][1])==FALSE){
        if(Set[[kk+1]][1]<mu.min){Set2[[kk+1]][1]=mu.min}
        if(Set[[kk+1]][length(Set[[kk+1]])]>mu.max){Set2[[kk+1]][length(Set[[kk+1]])]=mu.max}}
      delta=.4 #size of "fleck"
      delta2=.005 #size of space between "flecks"
      for(ii in 1:mm) {lines(Set2[[kk+1]][2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1])),2),lwd=3,col=co[kk==OLOC])
                       #browser()
                       lines(rep(Set2[[kk+1]][2*ii-1],2)+delta2,c(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))-delta,min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))+delta),lwd=3,col=co[kk==OLOC])
                       lines(rep(Set2[[kk+1]][2*ii-0],2)-delta2,c(min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))-delta,min(c[OLOC+1]-0.25*b[OLOC+1]^2/(1e-3+a[OLOC+1]))+delta),lwd=3,col=co[kk==OLOC])
                       ######boldbit
                       if(is.na(Set2[[kk+1]][2*ii-1])==FALSE){
                         x2<-seq(Set2[[kk+1]][2*ii-1],Set2[[kk+1]][2*ii-0],length=50)
                         lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[kk==OLOC],lwd=4)}
                       ######
      }
    }
  }
  legend("topleft", legend=paste(OLOC," "), col=co[OLOC%in%OLOCdash],pch=15,bg="white")
  ####PLOT END########
  ######PLOT timestep end##########
  OLOCdash2<-c()
  for(it in OLOC){
    if(is.na(Set[[it+1]][1])==FALSE){
      OLOCdash2<-c(OLOCdash2,it)}
  }
  for(kk in OLOCdash2){
    if(kk==OLOCdash2[1]){
      plot(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],type="l",col=co[kk==OLOCdash],ylim=F[taustar+1]+c(0,10),xlab=expression(mu),ylab="cost",xlim=c(mu.min,mu.max))
    }else{
      lines(x,a[kk+1]*x^2+b[kk+1]*x+c[kk+1],col=co[kk==OLOCdash])
    }
    if(length(Set[[kk+1]])>0){
      mm=length(Set[[kk+1]])/2
      Set2<-Set
      if(is.na(Set[[kk+1]][1])==FALSE){
        if(Set[[kk+1]][1]<mu.min){Set2[[kk+1]][1]=mu.min}
        if(Set[[kk+1]][length(Set[[kk+1]])]>mu.max){Set2[[kk+1]][length(Set[[kk+1]])]=mu.max}}
      delta=.4 #size of "fleck"
      delta2=.005 #size of space between "flecks"
      for(ii in 1:mm) {lines(Set2[[kk+1]][2*ii-c(1,0)]+c(delta2,-delta2),rep(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1])),2),lwd=3,col=co[kk==OLOCdash])
                       #browser()
                       lines(rep(Set2[[kk+1]][2*ii-1],2)+delta2,c(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))-delta,min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))+delta),lwd=3,col=co[kk==OLOCdash])
                       lines(rep(Set2[[kk+1]][2*ii-0],2)-delta2,c(min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))-delta,min(c[OLOCdash2+1]-0.25*b[OLOCdash2+1]^2/(1e-3+a[OLOCdash2+1]))+delta),lwd=3,col=co[kk==OLOCdash])
                       ######boldbit
                       if(is.na(Set2[[kk+1]][2*ii-1])==FALSE){
                         x2<-seq(Set2[[kk+1]][2*ii-1],Set2[[kk+1]][2*ii-0],length=50)
                         lines(x2,a[kk+1]*x2^2+b[kk+1]*x2+c[kk+1],col=co[kk==OLOCdash],lwd=4)}
                       ######
      }
    }
  }
  legend("topleft", legend=paste(OLOCdash2," "), col=co[OLOCdash%in%OLOCdash2],pch=15,bg="white")
  ####PLOT END########
    #print(taustar)
}
  }
  #return(c(F[n+1],CP[[n+1]]))
}


####Code to run the function####
#y<-SimChange(1000,c(1,2,3,0,5,4,6,5,1),1.5)
#plot(y[[2]],typ="l")
#abline(v=y[[1]],col="red")

#OPR<-OptimalPartitioningRig(y[[2]],negloglike.norm.mean2,log(1000));OPR
#abline(v=OPR,col="blue")

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



