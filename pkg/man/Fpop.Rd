\name{Fpop}
\alias{Fpop}
\title{Fpop}
\description{Function calling the fpop algorithm, use functional pruning and
optimal partionning to recover the best segmentation with respect
to the L2 loss with a per change-point penalty of lambda. More
precisely, this function computes the solution to argmin_m
sum_{i=1}^n (x_i-m_i)^2 + lambda * sum_{i=1}^{n-1} I(m_i !=
m_{i+1}), where the indicator function I counts the number of
changes in the mean vector m.}
\usage{Fpop(x, lambda, mini = min(x), maxi = max(x))}
\arguments{
  \item{x}{A vector of double : the signal to be segmented}
  \item{lambda}{Value of the penalty}
  \item{mini}{Min value for the mean parameter of the segment}
  \item{maxi}{Max value for the mean parameter of the segment}
}

\value{return a list with a vector t.est containing the position of the change-points}

\author{Guillem Rigaill, Toby Dylan Hocking}




\examples{
set.seed(1)
N <- 100
data.vec <- c(rnorm(N), rnorm(N, 2), rnorm(N))
fit <- Fpop(data.vec, N)
end.vec <- fit$t.est
change.vec <- end.vec[-length(end.vec)]
start.vec <- c(1, change.vec+1)
segs.list <- list()
for(seg.i in seq_along(start.vec)){
  start <- start.vec[seg.i]
  end <- end.vec[seg.i]
  segs.list[[seg.i]] <- data.frame(
    start, end, mean=mean(data.vec[start:end]))
}
segs <- do.call(rbind, segs.list)
plot(data.vec)
with(segs, segments(start-0.5, mean, end+0.5, mean, col="green"))
with(segs[-1,], abline(v=start-0.5, col="green", lty="dotted"))
}
