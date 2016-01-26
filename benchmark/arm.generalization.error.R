source("packages.R")

source("pick.best.index.R")

load("all.stats.RData")

estimate.test.error <- function
### Do leave-one-out cross-validation on chromosome arms.
(stats
### Named list with arrays errors, false.positive, false.negative,
### each of dim nparam x nprof x nfolds.
 ){
  stats <- stats[c("errors","false.positive","false.negative")]
  nparam <- dim(stats$errors)[1]
  nprof <- dim(stats$errors)[2]
  nfolds <- dim(stats$errors)[3]
  stat.names <- names(stats)
  local.loo <- matrix(NA,length(stats),nfolds)
  hybrid.loo <- matrix(NA,length(stats),nfolds)
  global.loo <- matrix(NA,length(stats),nfolds)
  rownames(global.loo) <- stat.names
  rownames(hybrid.loo) <- stat.names
  rownames(local.loo) <- stat.names
  train.err.mat <-
    matrix(NA,3,nfolds,dimnames=list(method=c("global","hybrid","local"),fold=NULL))
  for(fold in 1:nfolds){

    train.err <- rep(NA,nparam) ## global model
    for(j in 1:nparam){
      train.err[j] <- sum(stats$errors[j,,-fold],na.rm=TRUE)
    }
    ## save for hybrid approach
    global.train.err <-
      data.frame(train.err,param=1:length(train.err))
    global.picked <- pick.best.index(train.err)
    for(sn in stat.names){
      global.loo[sn,fold] <-
        mean(stats[[sn]][global.picked,,fold],na.rm=TRUE)
    }
    train.err.mat["global",fold] <- train.err[global.picked] ## for comparing train err

    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    hybrid.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## hybrid models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      is.min <- train.err == min(train.err)
      global.subset <- global.train.err[is.min,]
      hybrid.picked <-
        with(global.subset,param[which.min(train.err)])
      hybrid.train.errors[i] <- train.err[hybrid.picked]
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][hybrid.picked,i,fold]
      }
    }
    hybrid.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["hybrid",fold] <- sum(hybrid.train.errors)
    
    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    local.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## local models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      local.picked <- pick.best.index(train.err)
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][local.picked,i,fold]
      }
      local.train.errors[i] <- train.err[local.picked]
    }
    local.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["local",fold] <- sum(local.train.errors)
    
  }
  list(local=local.loo,
       hybrid=hybrid.loo,
       global=global.loo,
       train.err.mat=train.err.mat)
### Named list with elements local, hybrid, global, each a 3 x nfolds
### matrix.
}

## just calculate errors with respect to fit!
stat.dfs <- lapply(names(all.stats),function(algorithm){
  one <- all.stats[[algorithm]]
  global.fp <- apply(one$false.positive,1,sum,na.rm=TRUE)
  global.fn <- apply(one$false.negative,1,sum,na.rm=TRUE)
  num.normal <- sum(one$normal.anns)
  num.breakpoint <- sum(one$breakpoint.anns)
  num.total <- num.normal+num.breakpoint
  fpfn.df <- function(false.positive,false.negative,type,parameter){
    FP <- false.positive/num.normal * 100
    FN <- false.negative/num.breakpoint * 100
    error <- (false.positive+false.negative)/num.total * 100
    data.frame(error,FP,FN,
               false.positive=false.positive*100/num.total,
               false.negative=false.negative*100/num.total,
               algorithm,
               parameter,
               type,
               row.names=NULL)
  }
  global.df <- fpfn.df(global.fp,global.fn,"global",one$parameter)
  ## nparam x nprof
  local.err.mat <- apply(one$error,1:2,sum,na.rm=TRUE)
  local.pick <- apply(local.err.mat,2,pick.best.index)
  get.fp.or.fn <- function(var){
    sum(sapply(seq_along(local.pick),function(i){
      param <- local.pick[i]
      sum(one[[var]][param,i,],na.rm=TRUE)
    }))
  }
  local.fp <- get.fp.or.fn("false.positive")
  local.fn <- get.fp.or.fn("false.negative")
  local.df <- fpfn.df(local.fp,local.fn,"local",NA)
  rbind(local.df,global.df)
})
## train error of best model, local and global.
best.global <- do.call(rbind,lapply(stat.dfs,function(df){
  global <- subset(df,type == "global")
  i <- which.min(global$error)
  global[i,]
}))
best.local <- do.call(rbind,lapply(stat.dfs,subset,type=="local"))
getmat <- function(df){
  m <- df[,c("error","FP","FN")]
  rownames(m) <- df$algorithm
  for(j in 1:ncol(m)){
    m[,j] <- sprintf("%.1f",m[,j])
  }
  as.matrix(m)
}

results <- lapply(all.stats,estimate.test.error)
## compare local, hybrid, and global error rates
print(sapply(results,function(stat)sapply(stat,function(x)mean(x[1,]))))
arm.generalization.error <- laply(results,"[[","global")
names(dimnames(arm.generalization.error)) <- c("algorithm","statistic","chromosome")
dimnames(arm.generalization.error)[[1]] <- names(all.stats)

## average over all arms
global.local.stats <- laply(results,function(L)do.call(rbind,lapply(L,rowMeans)))
names(dimnames(global.local.stats)) <- c("algorithm","method","statistic")
dimnames(global.local.stats)[[1]] <- names(all.stats)
dimnames(global.local.stats)[[3]] <- c("errors","FP","FN")

## now figure out a good ordering
arm.error <- arm.generalization.error[,1,]
ord <- rownames(arm.error)[order(apply(arm.error,1,mean))]

## this will be attached to the bigger table later.
train.err.mat <- cbind(getmat(best.global),getmat(best.local))[ord,]


#### figure: training error curves.
## figure out which is the global optimum.
id.vars <- c("param.id","parameter","algorithm")
measure.vars <- c("error","false.positive","false.negative")
global.train <- do.call(rbind,lapply(stat.dfs,function(df){
  plot.cols <- names(df)[names(df)%in%c(id.vars,measure.vars)]
  only.global <- subset(df,type=="global")[,plot.cols]
  if(nrow(only.global)==1){
    data.frame() ## exclude algos with no tuning parameter
  }else{
    only.global$param.id <- 1:nrow(only.global)
    melt(only.global,id=id.vars)
  }
}))
global.train$algorithm <- factor(global.train$algorithm,ord)
## minimize each curve
daply(global.train,.(algorithm),function(df)range(df$parameter))
min.vals <- ddply(global.train,.(algorithm),function(df){
  error <- subset(df,variable=="error")
  i <- which.min(error$value)
  error[i,]
})
min.vals$label <- sprintf("%.1f",min.vals$value)

train.err.curves <- ggplot(global.train,aes(param.id,value))+
  geom_line(aes(group=variable,colour=variable,linetype=variable),
            lwd=2)+
  geom_segment(aes(yend=value),xend=1,data=min.vals)+
  geom_text(aes(label=label),x=1,data=min.vals,vjust=-0.5,hjust=0)+
  facet_grid(.~algorithm,scales="free_x",space="free")+
  scale_colour_manual(values=c(
                        false.positive="grey50",
                        false.negative="#E41A1C",
                        error="black"
                        ))+
  scale_linetype_manual(values=c(error=1,false.positive=2,false.negative=3))+
  theme_bw()+ # bioinformatics
  theme(panel.margin=grid::unit(0,"lines"))+
  ylab("percent incorrectly predicted annotations in training set")
algos <- unique(global.train$algo)
n.algos <- length(algos)
png("figure-train-error-curves.png",height=400,width=300*n.algos)
print(train.err.curves)
dev.off()

## times figure
png("figure-profile-size-time.png",height=400,width=300*n.algos)
data(neuroblastoma,package="neuroblastoma")
profile.list <- split(neuroblastoma$profiles,neuroblastoma$profiles$profile.id)
all.cids <- names(all.stats[[1]]$seconds)
probes <- sapply(profile.list[all.cids], nrow)
seconds <- do.call(rbind,lapply(names(all.stats),function(model){
    L <- all.stats[[model]]
      data.frame(seconds=L$seconds,model,probes)
  }))
fastest.first <- sort(sapply(all.stats,function(L)median(L$seconds)))

seconds$model <- factor(seconds$model,names(fastest.first))
median.df <- data.frame(model=factor(names(fastest.first)),seconds=fastest.first)

good.breaks <- c(min(probes),max(probes[probes!=max(probes)]),max(probes))
time.labels <- do.call(rbind,lapply(names(all.stats),function(model){
  seconds <- c("1 second"=1,"1 minute"=60,"1 hour"=60*60,"1 day"=60*60*24)
  data.frame(model=factor(model,names(fastest.first)),
               seconds,label=names(seconds),probes=20000)
}))
p <- ggplot()+
  geom_text(aes(probes,seconds,label=label),data=time.labels)+
  geom_point(aes(probes,seconds),data=seconds,pch=1)+
  geom_hline(aes(yintercept=seconds),data=median.df)+
  facet_grid(.~model)+
  ggtitle("Algorithms show varying speed and complexity")+
  theme(panel.margin=grid::unit(1,"lines"))+
  scale_x_log10("Number of probes on the profile",
                breaks=good.breaks,minor_breaks=NULL)+
  scale_y_log10("Seconds to process the profile")

print(p)
dev.off()

## old version (divide by total number of examples)
fpfn <- cbind(global.local.stats[,"global",],
              global.local.stats[,"local",])
## new version (divide by number of positive/neg examples)
num.normal <- do.call(rbind,lapply(all.stats,function(stat){
  colSums(stat$normal.anns)
}))
num.breakpoint <- do.call(rbind,lapply(all.stats,function(stat){
  colSums(stat$breakpoint.anns)
}))
num.total <- num.normal+num.breakpoint
fp.count <- arm.generalization.error[,"false.positive",]*num.total
stopifnot(all(abs(round(fp.count)-fp.count)<1e-6))
fp <- fp.count/num.normal
fn <- arm.generalization.error[,"false.negative",]*num.total/num.breakpoint
## do for each local and global
fpfn.list <- list()
for(training.method in c("global","local")){
  errors <- rep(NA,length(results))
  names(errors) <- names(results)
  FP <- rep(NA,length(results))
  names(FP) <- names(results)
  FN <- rep(NA,length(results))
  names(FN) <- names(results)
  for(algorithm in names(results)){
    m <- results[[algorithm]][[training.method]]
    FP[algorithm] <-
      mean(m["false.positive",]*num.total[algorithm,]/num.normal[algorithm,])
    FN[algorithm] <-
     mean(m["false.negative",]*num.total[algorithm,]/num.breakpoint[algorithm,])
    errors[algorithm] <- mean(m["errors",])
  }
  fpfn.list[[training.method]] <- cbind(errors,FP,FN)
}
fpfn <- with(fpfn.list,cbind(global,local))[ord,]*100
fpfn[1:length(fpfn)] <- sprintf("%3.1f",fpfn)

save(arm.generalization.error,
     num.normal,
     num.breakpoint,
     num.total,
     file="arm.generalization.error.RData")
