source("packages.R")

##system("sudo apt-get install libgsl0-dev")

## the job of this script is to make zzz.stats.RData, which caches the
## results of the model smoothing in a nice R data structure
smoothdir <- "smooth"

data(neuroblastoma)

source("smoothing-functions.R")
algos <- c(
  "smuce.penalized", "wbs.penalized",
  "multiBinSeg",
  "pelt.n",
  "cghseg.k",
  "fpop")

profiles.by.id <- with(neuroblastoma, {
  split(profiles, profiles$profile.id, drop=TRUE)
})
annotations.by.id <- with(neuroblastoma, {
  split(annotations, annotations$profile.id, drop=TRUE)
})
stopifnot(identical(names(profiles.by.id), names(annotations.by.id)))
all.cids <- names(profiles.by.id)
all.cids <- names(profiles.by.id)[1:100]
chrom.order <- as.character(c(1,2,3,4,11,17))

## Run smoothers for all algos.
for(clin.id.i in seq_along(all.cids)){
  clin.id <- all.cids[[clin.id.i]]
  cat(sprintf("%4d / %4d profile=%s\n",
              clin.id.i, length(profiles.by.id), clin.id))
  one <- profiles.by.id[[clin.id]]
  some.chroms <- subset(one, chromosome %in% chrom.order)
  these.labels <- annotations.by.id[[clin.id]]
  print(these.labels)
  for(a in algos){
    errors.csv.gz <- file.path(smoothdir, clin.id, a, "errors.csv.gz")
    if(!file.exists(errors.csv.gz)){
      update.smoothers <- smoothers[a]
      run.smoothers(one, these.labels, update.smoothers, db=smoothdir)
    }
  }
}

all.stats <- list()
## each all.stats array is nparam x nprofiles x nann
for(a in algos){
  print(a)
  f <- file.path(smoothdir,clin.id,a,"parameters.csv.gz")
  parameters <- tryCatch({
    scan(f,quiet=TRUE)
  },error=function(e){
    scan(f,quiet=TRUE,what="char")
  })
  param.names <- as.character(parameters)
  breakpoint.anns <-
    matrix(0,length(all.cids),length(chrom.order),
           dimnames=list(profile=all.cids,chromosome=chrom.order))
  normal.anns <- breakpoint.anns
  errors <-
    array(NA,
          list(length(param.names),length(all.cids),length(chrom.order)),
          list(param.names,all.cids,chrom.order))
  labels <- errors
  predictions <- errors
  false.positive <- errors
  false.negative <- errors
  for(cid in all.cids){
    errfile <- file.path(smoothdir,cid,a,"errors.csv.gz")
    e <- as.matrix(read.csv(errfile,header=FALSE))
    labfile <- file.path(smoothdir,cid,a,"breakpoint.labels.csv.gz")
    anns <-
      read.csv(labfile,header=FALSE,
               col.names=c("profile.id","chromosome","min","max","annotation"))
    ann.mat <- matrix(as.character(anns$annotation),
                      nrow=nrow(e),ncol=ncol(e),byrow=TRUE)
    colnames(e) <- anns$chromosome
    errors[,cid,colnames(e)] <- e
    for(chr in colnames(e)){
      normal.anns[cid,chr] <- if(is.na(e[1,chr]))0 else
        nrow(subset(anns,chromosome==chr & annotation=="normal"))
      breakpoint.anns[cid,chr] <- if(is.na(e[1,chr]))0 else
        nrow(subset(anns,chromosome==chr & annotation=="breakpoint"))
    }
    ## Careful: is.na(NA & TRUE) but !is.na(NA & FALSE)
    false.negative[,cid,colnames(e)] <-
      ifelse(is.na(e),NA,ifelse(e & ann.mat=="breakpoint",1L,0L))
    false.positive[,cid,colnames(e)] <-
      ifelse(is.na(e),NA,ifelse(e & ann.mat=="normal",1L,0L))
  }
  ## TODO: use "finished" file instead of checking for all these
  ## files.
  readsecs <- function(cid){
    secfile <- file.path(smoothdir,cid,a,"seconds.csv.gz")
    scan(secfile,quiet=TRUE)
  }
  seconds <- sapply(all.cids,readsecs)
  all.stats[[a]] <- list(errors=errors,
                         false.positive=false.positive,
                         false.negative=false.negative,
                         parameters=parameters,
                         seconds=seconds,
                         normal.anns=normal.anns,
                         breakpoint.anns=breakpoint.anns)
}

save(all.stats, file="all.stats.RData")
