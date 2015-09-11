works_with_R("3.2.2", xtable="1.7.4", ggplot2="1.0.1",
             tikzDevice="0.7.0")

options(tikzDocumentDeclaration="\\documentclass{article}\\usepackage{amssymb,amsmath}",
        tikzMetricsDictionary="tikzMetrics")

if(!file.exists("arm.generalization.error.RData")){
  download.file("http://cbio.ensmp.fr/~thocking/neuroblastoma/arm.generalization.error.RData", "arm.generalization.error.RData")
}
load("arm.generalization.error.RData")

if(!file.exists("all.stats.RData")){
  download.file("http://cbio.ensmp.fr/~thocking/neuroblastoma/zzz.stats.RData", "all.stats.RData")
}
load("all.stats.RData")

algo.colors <-
  c(pDPA="#1B9E77",
    pelt="#D95F02", pelt.default="#D95F02",
    fpop="#7570B3",
    binseg="#E7298A",
    ##dnacopy.default="#66A61E",
    SNIP="#E6AB02",
    wbs="#A6761D", wbs.default="#A6761D",
    smuce="#666666", smuce.default="#666666")

algo.sizes <-
  c(pDPA=1,
    pelt=1, pelt.default=1,
    fpop=1,
    binseg=1,
    ##dnacopy.default="#66A61E",
    ##SNIP="#E6AB02",
    wbs=2, wbs.default=2,
    smuce=2, smuce.default=2)

algos <-
  c(pelt="pelt.n",
    pDPA="cghseg.k",
    fpop="fpop",
    binseg="multiBinSeg",
    wbs="wbs.th.const",
    pelt.default="pelt.default",
    smuce="smuce.alpha",
    smuce.default="smuce.default",
    wbs.default="wbs.default")

percent.mat <- arm.generalization.error[algos, "errors", ] * 100
rownames(percent.mat) <- names(algos)

percent.df <-
  data.frame(mean=rowMeans(percent.mat),
             sd=apply(percent.mat, 1, sd))

percent.wide <-
  data.frame(percent.df,
             algorithm=factor(rownames(percent.df), rownames(percent.df)),
             min=apply(percent.mat, 1, min),
             max=apply(percent.mat, 1, max))

dot.algos <- rownames(percent.df)[percent.df$mean < 10]
dot.algos <- rownames(percent.df)
percent.by.algo <- list()
roc.by.algo <- list()
min.err.by.algo <- list()
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
  print(fp.possible + fn.possible)
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
               row.names=NULL)
  if(!algorithm %in% c("binseg", "pDPA", "fpop")){
    roc.by.algo[[algorithm]] <- algo.df
  }
  min.param <- algo.df[which.min(errors), ]
  if(algorithm %in% c("wbs", "smuce")){
    min.err.by.algo[[algorithm]] <- min.param
  }
  if(algorithm == "pelt"){
    MLseg <- min.param
  }
}
min.err <- do.call(rbind, min.err.by.algo)
all.roc <- do.call(rbind, roc.by.algo)
roc.curves <- subset(all.roc, !grepl("default", algorithm))
default.dots <- subset(all.roc, grepl("default", algorithm))
roc.labels <-
  rbind(data.frame(default.dots, param="default"),
        data.frame(min.err, param="min error"))
roc.dots <-
  rbind(roc.labels,
        data.frame(MLseg, param="min error"))
percent.tall <- do.call(rbind, percent.by.algo)
rownames(roc.labels) <- roc.labels$algorithm
roc.labels["wbs", c("fpr.percent", "tpr.percent")] <- c(5, 88)
roc.labels["pelt.default", c("fpr.percent")] <- 5
roc.labels["smuce", c("tpr.percent")] <- 95
roc.labels["smuce.default", c("fpr.percent", "tpr.percent")] <- c(73, 95)
roc.labels["wbs.default", c("fpr.percent", "tpr.percent")] <- c(90, 95)

tsize <- 3
ggroc <- 
ggplot()+
  scale_x_continuous("False positive rate (percent)", limits=c(-15, 100))+
  scale_y_continuous("True positive rate (percent)",
                     breaks=seq(0, 100, by=25),
                     limits=c(0, 110))+
  scale_color_manual(values=algo.colors)+
  scale_size_manual(values=algo.sizes)+
  scale_fill_manual("parameter",
                    values=c(default="black", "min error"="white"))+
  geom_path(aes(fpr.percent, tpr.percent,
                color=algorithm, size=algorithm),
                data=roc.curves)+
  geom_point(aes(fpr.percent, tpr.percent,
                 fill=param,
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
  geom_text(aes(-15, 110, label=sprintf("fpop=pelt=pDPA=binseg\nFP=%d\nFN=%d",
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
  scale_x_continuous("percent incorrect labels",
                     limits=c(-90, max(percent.tall$percent)),
                     breaks=seq(0, 80, by=20))+
  geom_text(aes(-5, algorithm,
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

tikz("figure-neuroblastoma-test-error.tex", w=3, h=2)
print(p)
dev.off()

xt <- xtable(percent.df[nrow(percent.df):1, ], digits=2)
print(xt, file="table-neuroblastoma-test-error.tex", floating=FALSE)
