works_with_R("3.2.2", xtable="1.7.4", ggplot2="1.0.1",
             tikzDevice="0.7.0")

options(tikzDocumentDeclaration="\\documentclass{article}\\usepackage{amssymb,amsmath}",
        tikzMetricsDictionary="tikzMetrics")

if(!file.exists("arm.generalization.error.RData")){
  download.file("http://cbio.ensmp.fr/~thocking/neuroblastoma/arm.generalization.error.RData", "arm.generalization.error.RData")
}
load("arm.generalization.error.RData")

if(!file.exists("all.stats.RData")){
  download.file("http://cbio.ensmp.fr/~thocking/neuroblastoma/all.stats.RData", "all.stats.RData")
}
load("all.stats.RData")

source("algo.colors.R")

algos <-
  c(BinSeg="multiBinSeg",
    PELT="pelt.n",
    pDPA="cghseg.k",
    FPOP="fpop",
    PELT.default="pelt.default",
    SMUCE="smuce.penalized",
    SMUCE.default="smuce.default",
    WBS.default="wbs.default",
    WBS="wbs.penalized")

percent.mat <- arm.generalization.error[algos, "errors", ] * 100
rownames(percent.mat) <- names(algos)

percent.df <-
  data.frame(mean=rowMeans(percent.mat),
             sd=apply(percent.mat, 1, sd))

dot.algos <- grep("default", rownames(percent.df), value=TRUE,invert=TRUE)
percent.wide <-
  data.frame(percent.df[dot.algos, ],
             algorithm=factor(dot.algos, dot.algos),
             min=apply(percent.mat[dot.algos, ], 1, min),
             max=apply(percent.mat[dot.algos, ], 1, max))
percent.by.algo <- list()
roc.by.algo <- list()
min.err.by.algo <- list()
breakpoint.tab.list <- list()
normal.tab.list <- list()
for(algorithm in dot.algos){
  percent.by.algo[[algorithm]] <-
    data.frame(algorithm=factor(algorithm, rownames(percent.df)),
               percent=percent.mat[algorithm, ])
  orig.algo <- algos[[algorithm]]
  algo.stats <- all.stats[[orig.algo]]
  fp <- apply(algo.stats$false.positive, 1, sum, na.rm=TRUE)
  fn <- apply(algo.stats$false.negative, 1, sum, na.rm=TRUE)
  errors <- fp + fn
  fp.possible <- sum(algo.stats$normal.anns)
  fn.possible <- sum(algo.stats$breakpoint.anns)
  breakpoint.tab.list[[algorithm]] <-
    table(algo.stats$breakpoint.anns, useNA="always")
  normal.tab.list[[algorithm]] <-
    table(algo.stats$normal.anns, useNA="always")
  total.labels <- fp.possible + fn.possible
  fpr.percent <- fp/fp.possible * 100
  fnr <- fn/fn.possible
  tpr.percent <- (1-fnr) * 100
  algo.df <- 
    data.frame(algorithm,
               parameter.i=seq_along(algo.stats$parameters),
               parameter=paste(algo.stats$parameters),
               fp,
               fn,
               fpr.percent,
               tpr.percent,
               total.labels,
               row.names=NULL)
  if(!algorithm %in% c("BinSeg", "pDPA", "PELT")){
    roc.by.algo[[algorithm]] <- algo.df
  }
  min.param <- algo.df[which.min(errors), ]
  if(algorithm %in% c("WBS", "SMUCE")){
    min.err.by.algo[[algorithm]] <- min.param
  }
  if(algorithm == "FPOP"){
    MLseg <- min.param
  }
}
(breakpoint.tab <- do.call(rbind, breakpoint.tab.list))
(normal.tab <- do.call(rbind, normal.tab.list))

with(all.stats, wbs.th.const$normal.anns - wbs.default$normal.anns)
head(all.stats$wbs.th.const$normal.anns)
head(all.stats$wbs.default$normal.anns)
with(all.stats, {
  rbind(wbs.th.const.normal=wbs.th.const$normal.anns["6",],
        wbs.default.normal=wbs.default$normal.anns["6",],
        wbs.th.const.breakpoint=wbs.th.const$breakpoint.anns["6",],
        wbs.default.breakpoint=wbs.default$breakpoint.anns["6",])
})

min.err <- do.call(rbind, min.err.by.algo)
all.roc <- do.call(rbind, roc.by.algo)
with(all.roc, stopifnot(total.labels[1] == total.labels))
roc.curves <- subset(all.roc, !grepl("default", algorithm))
##default.dots <- subset(all.roc, grepl("default", algorithm))
roc.labels <-
  rbind(##data.frame(default.dots, param="default"),
        data.frame(min.err, param="min error"))
roc.dots <-
  rbind(roc.labels,
        data.frame(MLseg, param="min error"))
percent.tall <- do.call(rbind, percent.by.algo)
rownames(roc.labels) <- roc.labels$algorithm
roc.labels["WBS", c("fpr.percent", "tpr.percent")] <- c(3, 88)
##roc.labels["PELT.default", c("fpr.percent")] <- 5
roc.labels["SMUCE", c("fpr.percent")] <- -6
##roc.labels["SMUCE.default", c("fpr.percent", "tpr.percent")] <- c(73, 95)
##roc.labels["WBS.default", c("fpr.percent", "tpr.percent")] <- c(90, 95)

tsize <- 3
ggroc <- 
ggplot()+
  scale_x_continuous("False positive rate (percent)",
                     breaks=seq(0, 100, by=5))+
  scale_y_continuous("True positive rate (percent)",
                     breaks=seq(0, 100, by=5))+
  coord_equal(ylim=c(75, 102), xlim=c(-7, 25))+
  scale_color_manual(values=algo.colors)+
  scale_fill_manual("parameter",
                    values=c(default="black", "min error"="white"))+
  geom_path(aes(fpr.percent, tpr.percent,
                color=algorithm),
              data=roc.curves)+
  geom_point(aes(fpr.percent, tpr.percent,
                 ##fill=param,
                 color=algorithm),
             data=roc.dots,
             pch=21, size=5)+
  geom_text(aes(fpr.percent, tpr.percent, color=algorithm,
                label=sprintf("%s\nFP=%d\nFN=%d",
                  sub("[.]", "\n", algorithm), fp, fn)),
            data=roc.labels,
            hjust=0,
            vjust=1,
            size=tsize)+
  guides(color="none", size="none")+
  geom_text(aes(-5, 100, label=sprintf("fpop=pelt=pDPA=binseg\nFP=%d\nFN=%d",
                           fp, fn),
                color=algorithm),
            size=tsize,
            hjust=0,
            vjust=1,
            data=MLseg)+
  theme_bw()+
  theme(legend.position=c(0.7, 0.3))

pdf("figure-neuroblastoma-roc.pdf", 5, 3)
print(ggroc)
dev.off()

is.right <- rownames(percent.wide) != "wbs.default"
right.labels <- percent.wide[is.right, ]
left.labels <- percent.wide[!is.right, ]

p <- ggplot()+
  theme_bw()+
  guides(color="none")+
  scale_x_continuous("test error, percent incorrect labels",
                     limits=c(-5, max(percent.tall$percent)),
                     breaks=seq(0, 20, by=5))+
  geom_text(aes(0, algorithm,
                color=algorithm,
                label=sprintf("$%.2f\\%% \\pm %05.2f$", mean, sd)),
            data=percent.wide,
            hjust=1,
            size=3)+
  ## geom_text(aes(max, algorithm,
  ##               label=sprintf("\\ $%.2f\\%% \\pm %.2f$", mean, sd)),
  ##           data=right.labels,
  ##           size=3,
  ##           hjust=0)+
  ## geom_text(aes(min, algorithm,
  ##               label=sprintf("$%.2f\\%% \\pm %.2f$\\ \\ ", mean, sd)),
  ##           data=left.labels,
  ##           size=3,
  ##           hjust=1)+
  geom_point(aes(percent, algorithm, color=algorithm),
             pch=1,
             data=percent.tall)+
  scale_color_manual(values=algo.colors)

tikz("figure-neuroblastoma-test-error.tex", w=3.2, h=1.8)
print(p)
dev.off()

t.test(percent.mat["FPOP", ], percent.mat["SMUCE", ],
       alternative="less", paired=TRUE)
t.test(percent.mat["FPOP", ], percent.mat["WBS", ],
       alternative="less", paired=TRUE)

xt <- xtable(percent.df[nrow(percent.df):1, ], digits=2)
print(xt, file="table-neuroblastoma-test-error.tex", floating=FALSE)
