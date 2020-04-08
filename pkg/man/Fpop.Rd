\name{Fpop}
\alias{Fpop}
\title{Fpop}
\description{Function calling the fpop algorithm, use functional pruning and
optimal partionning to recover the best segmentation with respect
to the L2 loss with a per change-point penalty of \code{lambda}. More
precisely, this function computes the solution to argmin_m
sum_{i=1}^n (x_i-m_i)^2 + \code{lambda} * sum_{i=1}^{n-1} [1+I(m_i !=
m_{i+1})], where the indicator function I counts the number of
changes in the mean vector m, so the total penalty term is
equivalent to the number of segments.}
\usage{Fpop(x, lambda, mini = min(x), 
    maxi = max(x))}
\arguments{
  \item{x}{A vector of double : the signal to be segmented}
  \item{lambda}{Value of the penalty}
  \item{mini}{Min value for the mean parameter of the segment}
  \item{maxi}{Max value for the mean parameter of the segment}
}

\value{Named list with the following elements: input data (signal, n,
\code{lambda}, min, max), path (best previous segment end up to each data
point), cost (optimal penalized cost up to each data point, where
loss is the square loss without the constant term and the penalty
is \code{lambda} times the number of segments), t.est (vector of overall
optimal segment ends), K (optimal number of segments), J.est
(total un-penalized cost of optimal model, aka total square
loss). To see how cost relates to J.est, see definition of J.est
in the R source code for this function. (last cost = J.est -
sum.of.squares + K*lambda)}

\author{Guillem Rigaill <guillem.rigaill@inra.fr> [aut, cre], Toby Hocking [aut], Robert Maidstone [aut], Michel Koskas [ctb], Paul Fearnhead [aut]}




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
  seg.data <- data.vec[start:end]
  seg.mean <- mean(seg.data)
  segs.list[[seg.i]] <- data.frame(
    start, end,
    mean=seg.mean,
    seg.cost=sum((seg.data-seg.mean)^2))
}
segs <- do.call(rbind, segs.list)
plot(data.vec)
with(segs, segments(start-0.5, mean, end+0.5, mean, col="green"))
with(segs[-1,], abline(v=start-0.5, col="green", lty="dotted"))
}
