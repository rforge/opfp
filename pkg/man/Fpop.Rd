\name{Fpop}
\alias{Fpop}
\title{Fpop}
\description{Function calling the fpop algorithm, use functional pruning and optimal partionning. It is implemented for the L2-loss functon}
\usage{Fpop(x, lambda, mini = min(x), maxi = max(x))}
\arguments{
  \item{x}{A vector of double : the signal to be segmented}
  \item{lambda}{Value of the penalty}
  \item{mini}{Min value for the mean parameter of the segment}
  \item{maxi}{Max value for the mean parameter of the segment}
}

\value{return a list with a vector containing the position of the change-points}

\author{Guillem Rigaill}





