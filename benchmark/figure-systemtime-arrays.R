works_with_R("3.1.0", dplyr="0.2", directlabels="2014.6.13", ggplot2="1.0",
             reshape2="1.2.2")

load("systemtime.arrays.RData")

refs <- data.frame(unit=c("1 second","1 minute"),
                   seconds=c(1, 60))
small.refs <- refs[1,]
abbrevs <- c(pelt.SIC="pelt", fpop.SIC="fpop",
             cghseg.52="pDPA", multiBinSeg.52="binseg")
timings <- systemtime.arrays %.%
  filter(!grepl("dnacopy", algorithm)) %.%
  mutate(models=ifelse(algorithm %in%
           c("fpop.SIC", "pelt.SIC", "dnacopy.default"),
           "one", "several"),
         abbrev=abbrevs[as.character(algorithm)])
counts <- table(timings$abbrev)
tit <-
  paste(length(counts), "algorithms on",
        counts[1],
        "tumor chromosome segmentation problems (system.time)")
algo.colors <-
  c(pDPA="#1B9E77", pelt="#D95F02", fpop="#7570B3", binseg="#E7298A",
    "#66A61E", "#E6AB02",  "#A6761D", "#666666")

wide <- dcast(timings, pid.chr ~ algorithm, value.var="seconds")
faster.counts <- wide %.%
  summarise(fpop.faster=sum(pelt.SIC > fpop.SIC),
            same.speed=sum(pelt.SIC == fpop.SIC),
            pelt.faster=sum(pelt.SIC < fpop.SIC))
faster.labels <-
  data.frame(x=c(-0.5, -1.5, -1.5), y=c(-1, -1, -2),
             label=c(sprintf("fpop faster\n%d problems",
               faster.counts$fpop.faster),
               sprintf("pelt faster\n%d problems",
                       faster.counts$pelt.faster),
               sprintf("fpop=pelt\n%d problems",
                       faster.counts$same.speed)))
ref.color <- "red"
segs <-
  data.frame(x=c(-1.5,-1.75),
             y=c(-1.15,-2.15),
             xend=c(-Inf, -2.35),
             yend=c(-2.35, -2.38))
pid.chr.pelt.fast <- as.character(filter(wide, pelt.SIC < 0.01)$pid.chr)
save(pid.chr.pelt.fast, file="pid.chr.pelt.fast.RData")

inaccurate <-
  data.frame(x=0.0115, y=0.005,
             label=paste("system.time inaccurate for",
               length(pid.chr.pelt.fast), "problems,",
               "     see microbenchmark figure"))
scatter <-
  ggplot()+
  geom_rect(aes(xmin=0, ymin=0,
                xmax=0.01, ymax=0.01),
            fill="violet", alpha=1/2)+
  geom_segment(aes(10^x, 10^y, xend=10^xend, yend=10^yend), data=segs,
               arrow=grid::arrow(type="closed", length=grid::unit(0.1, "in")))+
  geom_hline(aes(yintercept=seconds),
             data=small.refs, color=ref.color)+
  geom_vline(aes(xintercept=seconds),
             data=refs, color=ref.color)+
  geom_text(aes(seconds, 10^-0.5, label=unit),
             data=refs, color=ref.color, angle=90, vjust=-0.5)+
  geom_text(aes(1e-2, seconds, label=unit),
            data=small.refs, vjust=1.5, hjust=0, color=ref.color)+
  geom_text(aes(10^x, 10^y, label=label), data=faster.labels)+
  geom_point(aes(pelt.SIC, fpop.SIC), data=wide)+
  coord_equal(xlim=10^c(-2.5, 2), ylim=10^c(-2.5, 0.1))+
  geom_abline()+
  ggtitle("4467 tumor chromosome segmentation problems (system.time)")+
  scale_x_log10("seconds to compute 1 segmentation using pelt",
                minor_breaks=NULL,
                breaks=10^seq(-2, 1, by=1))+
  scale_y_log10("seconds to compute 1 segmentation using fpop",
                minor_breaks=NULL,
                breaks=10^seq(-2, 0, by=1))+
  theme_grey()+
  geom_text(aes(x, y, label=label),
            data=inaccurate, color="violet", hjust=0)

pdf("figure-systemtime-arrays-fpop-pelt.pdf", w=10)
print(scatter)
dev.off()

binseg.counts <- wide %.%
  summarise(fpop.faster=sum(multiBinSeg.52 > fpop.SIC),
            same.speed=sum(multiBinSeg.52 == fpop.SIC),
            binseg.faster=sum(multiBinSeg.52 < fpop.SIC))
binseg.labels <- with(binseg.counts, {
  data.frame(x=c(-1.75, -1.75, -1.75), y=c(-2.45, -0.75, -1.9),
             angle=c(0, 0, 45),
             label=c(sprintf("fpop faster\n%d problems",
               fpop.faster),
               sprintf("binary segmentation faster\n%d problems",
                       binseg.faster),
               sprintf("binary segmentation=fpop\n%d problems",
                       same.speed)))
})
binseg.inacc <- wide %.%
  filter(fpop.SIC < 0.01,
         multiBinSeg.52 < 0.01)
scatter2 <-
  ggplot()+
  geom_rect(aes(xmin=0, ymin=0,
                xmax=0.01, ymax=0.01),
            fill="violet", alpha=1/2)+
  geom_text(aes(0.0115, 0.006,
             label=paste("system.time inaccurate for",
               nrow(binseg.inacc), "problems")),
            color="violet", hjust=0)+
  geom_hline(aes(yintercept=seconds),
             data=small.refs, color=ref.color)+
  geom_vline(aes(xintercept=seconds),
             data=small.refs, color=ref.color)+
  geom_text(aes(seconds, 10^-0.5, label=unit),
             data=small.refs, color=ref.color, angle=90, vjust=-0.5)+
  geom_text(aes(1e-2, seconds, label=unit),
            data=small.refs, vjust=1.5, hjust=0, color=ref.color)+
  geom_text(aes(10^x, 10^y, label=label, angle=angle), data=binseg.labels)+
  geom_point(aes(multiBinSeg.52, fpop.SIC), data=wide)+
  coord_equal(xlim=10^c(-2.6, 0.1), ylim=10^c(-2.6, 0.1))+
  geom_abline()+
  ggtitle("4467 tumor chromosome segmentation problems (system.time)")+
  scale_x_log10("seconds to compute 1 segmentation using binary segmentation",
                minor_breaks=NULL)+
  scale_y_log10("seconds to compute 1 segmentation using fpop",
                minor_breaks=NULL)+                
  theme_grey()

pdf("figure-systemtime-arrays-fpop-multiBinSeg.pdf")
print(scatter2)
dev.off()

with.leg <-
  ggplot()+
  geom_hline(aes(yintercept=seconds), data=refs)+
  geom_text(aes(10^5.2, seconds, label=unit),
            data=refs, vjust=-0.5, hjust=0)+
  geom_point(aes(probes, seconds, color=abbrev,
                 shape=models),
             data=timings)+
  scale_shape_manual(values=c(one=1, several=19))+
  scale_color_manual(values=algo.colors)+
  theme_grey()+
  scale_x_log10("number of data points to segment",
                limits=10^c(1.3, 5.7),
                minor_breaks=NULL,
                breaks=c(10^seq(2, 4, by=1), range(timings$probes)))+
  scale_y_log10("seconds",labels=function(x)sprintf("%.2f", x),
                minor_breaks=NULL,
                breaks=10^seq(-2, 4, by=1))+
  ggtitle(tit)
print(with.leg)

pos.method <-
  list("last.points",
       "calc.boxes",
       dl.trans(h=h*3/2, x=x+0.1),
       "calc.borders",
       ##"draw.rects",
       qp.labels("y","bottom","top",make.tiebreaker("x","y"),ylimits))

dl <- direct.label(with.leg, pos.method)

pdf("figure-systemtime-arrays.pdf", w=10)
print(dl)
dev.off()

