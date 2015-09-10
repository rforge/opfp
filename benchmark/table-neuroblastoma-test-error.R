works_with_R("3.2.2", xtable="1.7.4", ggplot2="1.0.1",
             tikzDevice="0.7.0")

options(tikzDocumentDeclaration="\\documentclass{article}\\usepackage{amssymb,amsmath}",
        tikzMetricsDictionary="tikzMetrics")

if(!file.exists("arm.generalization.error.RData")){
  download.file("http://cbio.ensmp.fr/~thocking/neuroblastoma/arm.generalization.error.RData", "arm.generalization.error.RData")
}

load("arm.generalization.error.RData")

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
for(algorithm in dot.algos){
  percent.by.algo[[algorithm]] <-
    data.frame(algorithm=factor(algorithm, rownames(percent.df)),
               percent=percent.mat[algorithm, ])
}
percent.tall <- do.call(rbind, percent.by.algo)

is.right <- rownames(percent.wide) != "wbs.default"
right.labels <- percent.wide[is.right, ]
left.labels <- percent.wide[!is.right, ]

p <- ggplot()+
  theme_bw()+
  scale_x_continuous("percent incorrect labels",
                     limits=c(-90, max(percent.tall$percent)),
                     breaks=seq(0, 80, by=20))+
  geom_text(aes(-5, algorithm,
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
  geom_point(aes(percent, algorithm),
             pch=1,
             data=percent.tall)

tikz("figure-neuroblastoma-test-error.tex", w=3, h=2)
print(p)
dev.off()

xt <- xtable(percent.df[nrow(percent.df):1, ], digits=2)
print(xt, file="table-neuroblastoma-test-error.tex", floating=FALSE)
