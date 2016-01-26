source("packages.R")

data(neuroblastoma)
profile.list <- with(neuroblastoma, split(profiles, profiles$profile.id))

seconds.file.vec <- Sys.glob("smooth/*/*/seconds.csv.gz")
timings.list <- list()
for(seconds.file.i in seq_along(seconds.file.vec)){
  seconds.file <- seconds.file.vec[[seconds.file.i]]
  algo.dir <- dirname(seconds.file)
  algorithm <- basename(algo.dir)
  sample.dir <- dirname(algo.dir)
  sample.id <- basename(sample.dir)
  probes <- nrow(profile.list[[sample.id]])
  con <- gzfile(seconds.file, "r")
  seconds <- scan(con)
  close(con)
  timings.list[[seconds.file]] <-
    data.frame(sample.id, algorithm, probes, seconds)
}
timings <- do.call(rbind, timings.list)

all.profiles <- 
  data.frame(sample.id=names(profile.list),
             probes=sapply(profile.list, nrow))
subset(all.profiles, ! sample.id %in% timings$sample.id)

ggplot()+
  geom_point(aes(log10(probes), 0),
             data=all.profiles)+
  geom_text(aes(log10(probes), log10(seconds),
                color=algorithm, label=sample.id),
            data=timings)

ggplot()+
  geom_text(aes(log10(probes), log10(seconds),
                color=algorithm, label=sample.id),
            data=subset(timings, probes < 1e4))
