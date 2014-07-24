works_with_R("3.1.1", dplyr="0.2", directlabels="2014.6.13", ggplot2="1.0")

load("systemtime.simulationLarge.RData")

refs <- data.frame(unit=c( "10 seconds"),
                   seconds=c( 10))
timings <- systemtime.simulation %.%
  mutate(models=ifelse(algorithm %in%
           c("fpop", "pelt"),
           "one", "several"))
counts <- table(timings$algorithm)
tit <-
  paste(length(counts), "algorithms on",
        counts[1],
        "simulated segmentation problems (system.time)")
with.leg <-
  ggplot()+
  geom_hline(aes(yintercept=seconds), data=refs)+
  geom_text(aes(10^4, seconds, label=unit),
            data=refs, vjust=-0.5, hjust=1)+
  geom_point(aes(Ktrue, seconds, color=algorithm,
                 shape=models),
             data=timings)+
  scale_shape_manual(values=c(one=1, several=19))+
  scale_color_brewer(palette="Dark2")+
  theme_grey()+
  scale_x_log10("number of true of changes",
                breaks=10^seq(0, 5, by=1),
                limits=c(0.5, 2e4),
                minor_breaks=NULL)+
  scale_y_log10("seconds",
                minor_breaks=NULL,
                breaks=10^seq(-3,3, by=1))+
  ggtitle(tit)

pos.method <-
 list("first.points",
      "calc.boxes",
      ##dl.trans(x=max(x)),
      dl.trans(h=h*3/2, x=x-0.1),
      "calc.borders",
      ##"draw.rects",
      qp.labels("y","bottom","top",make.tiebreaker("x","y"),ylimits))

dl <- direct.label(with.leg, pos.method)

pdf("figure-systemtime-simulationLarge.pdf")
print(dl)
dev.off()


