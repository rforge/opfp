\name{Fpop}
\alias{Fpop}
\title{Fpop}
\description{Function calling the fpop algorithm, use functional pruning and
optimal partionning to recover the best segmentation with respect
to the L2 loss with a per change-point penalty of lambda. More
precisely, this function computes the solution to argmin_m
sum_{i=1}^n (x_i-m_i)^2 + lambda * sum_{i=1}^{n-1} I(m_i \neq
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

\author{Guillem Rigaill}





