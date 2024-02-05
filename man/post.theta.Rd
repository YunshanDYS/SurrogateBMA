\name{post.theta}
\alias{post.theta}
\title{
Generates posterior samples of the parameters
}
\description{
Provides a function to generate MCMC samples of p(theta_m | M = m, Data).
}
\usage{
post.theta(Y, S, A, m, nmc = 500, prior.para = NULL)
}
\arguments{
  \item{Y}{
}
  \item{S}{
}
  \item{A}{
}
  \item{m}{
}
  \item{nmc}{
}
  \item{prior.para}{
}
}

\value{
a numeric matrix; each row contains a MC sample from the posterior of the parameters given a certain model.
}
\references{
Duan, Y. and Parast, L., 2023. Flexible evaluation of surrogate markers with Bayesian model averaging. Statistics in Medicine.
}
\author{
Yunshan Duan
}

\keyword{internal}
