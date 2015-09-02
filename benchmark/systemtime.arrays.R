works_with_R("3.2.2",
             wbs="1.1",
             microbenchmark="1.4.2",
             SegAnnot="1.2",
             fpop="2014.7.16",
             Segmentor3IsBack="1.8",
             changepoint="1.1.2",
             DNAcopy="1.36.0",
             cghseg="1.0.2")

### Run several algorithms to perform a speed comparison:
### PrunedDP-cghseg/PrunedDP-segmentor/fpop/pelt/dnacopy/binary segmentation.
if(!file.exists("signal.list.RData")){
  u <- "http://cbio.ensmp.fr/~thocking/data/signal.list.RData"
  download.file(u, "signal.list.RData")
}
load("signal.list.RData")

run.wbs <- function(pro, th.const.vals){
  require(wbs)
  each.chrom(pro, function(chr){
    Y <- chr$logratio
    mean.mat <- matrix(NA, length(th.const.vals), length(Y))
    if(nrow(chr) >= 2){
      set.seed(1)
      w <- wbs(Y)
      list.or.na <- changepoints(w, th.const=th.const.vals)$cpt.th
      if(is.list(list.or.na)){
        for(param.i in seq_along(list.or.na)){
          changes <- sort(list.or.na[[param.i]])
          ends <- c(changes, length(Y))
          starts <- c(1, changes+1)
          for(seg.i in seq_along(starts)){
            start <- starts[[seg.i]]
            end <- ends[[seg.i]]
            mean.mat[param.i, start:end] <- mean(Y[start:end])
          }
        }
      }
    }
    mean.mat
  })
}

seg.funs <-
  list(wbs=function(one.chrom){
    wbs(one.chrom$logratio)
  },cghseg.52=function(one.chrom){
    kmax <- min(52, nrow(one.chrom))
    with(one.chrom, run.cghseg(logratio, position, kmax))
  },
  ##      dnacopy.sd=function(one.chrom){
  ##   probes <- nrow(one.chrom)
  ##   cna <- with(one.chrom, {
  ##     CNA(logratio, rep(1, probes), maploc=position)
  ##   })
  ##   sdvals <- sort(c(seq(20, 4, l = 20), 2^seq(2.5, -5, l = 10)[-1], 3, 0))
  ##   for(param.i in seq_along(sdvals)){
  ##     undo.SD <- sdvals[[param.i]]
  ##     segment(cna, undo.SD=undo.SD, undo.splits="sdundo",
  ##             verbose=0)
  ##   }
  ## },
       dnacopy.default=function(one.chrom){
    probes <- nrow(one.chrom)
    cna <- with(one.chrom, {
      CNA(logratio, rep(1, probes), maploc=position)
    })
    segment(cna, verbose=0)
  }, pelt.SIC=function(one.chrom){
    cpt.mean(one.chrom$logratio, method="PELT", penalty="Manual",
             pen.value=log(nrow(one.chrom)))
  }, fpop.SIC=function(one.chrom){
    Fpop(one.chrom$logratio, log(nrow(one.chrom)))
  },
## Segmentor gave a segfault:
## 1627 / 4467 44.19
## 1628 / 4467 57.19
## 1629 / 4467 64.19
## 1630 / 4467 67.19
## 1631 / 4467 123.19
## 1632 / 4467 145.19
## 1633 / 4467 162.19
## 1634 / 4467 173.19

##  *** caught segfault ***
## address (nil), cause 'unknown'

## Traceback:
##  1: .C("SegmentNormal", Size = as.integer(n2), KMax = as.integer(Kmax),     Data = as.double(dat), DataComp = as.integer(datasize), Breakpoints = as.integer(breaks),     Parameters = as.double(parameters), Likelihood = as.double(likelihood),     PACKAGE = "Segmentor3IsBack")
##  2: Segmentor.default(one.chrom$logratio, model = 2, Kmax = kmax)
##  3: Segmentor(one.chrom$logratio, model = 2, Kmax = kmax)
##  4: fun(one.chrom)
##  5: system.time({    fun(one.chrom)})
## aborting ...
## Segmentation fault (core dumped)
## make: *** [timings.RData] Error 139
## thocking@silene:~/breakpoints/webapp/applications-n
       
  ##      Segmentor3IsBack.52=function(one.chrom){
  ##   kmax <- min(52, nrow(one.chrom))
  ##   Segmentor(one.chrom$logratio, model=2, Kmax=kmax)
  ## },
       multiBinSeg.52=function(one.chrom){
    max.breaks <- min(51, nrow(one.chrom)-1) ## BREAKPOINTS not segments.
    multiBinSeg(one.chrom$logratio, max.breaks)
  })

## TODO: Figure with 2 panels. Several segmentations: multiBinSeg,
## Segmentor, cghseg, dnacopy.sd. One segmentation: fpop, pelt,
## multiBinSeg, dnacopy.default.

## TODO: scatterplot pelt vs fpop (coord_equal). Expect pelt > fpop.
    
systemtime.arrays.list <- list()
for(pid.chr.i in seq_along(signal.list)){
  pid.chr <- names(signal.list)[pid.chr.i]
  cat(sprintf("%4d / %4d %s\n", pid.chr.i, length(signal.list), pid.chr))
  one.chrom <- signal.list[[pid.chr]]
  m.args <- list(times=1)
  probes <- nrow(one.chrom)
  for(algorithm in names(seg.funs)){
    m.args[[algorithm]] <-
      substitute(fun(one.chrom), list(fun=seg.funs[[algorithm]]))
  }
  times <- do.call(microbenchmark, m.args)
  seconds <- times$time/1e9
  systemtime.arrays.list[[pid.chr]] <- 
    data.frame(algorithm=times$expr, pid.chr, probes, seconds)
}
systemtime.arrays <- do.call(rbind, systemtime.arrays.list)

save(systemtime.arrays, file="systemtime.arrays.RData")
