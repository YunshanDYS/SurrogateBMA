\name{post.model}
\alias{post.model}
\title{
Calculates the posterior probability of the candidate models
}
\description{
Gives the posterior probability of each candidate model being true.
}
\usage{
post.model(Y, S, A, prior.para = NULL)
}
\arguments{
  \item{Y}{
numeric vector; primary outcome, assumed to be continuous.
}
  \item{S}{
numeric vector; surrogate marker, assumed to be continuous.
}
  \item{A}{
numeric vector; treatment arm, assumed to be binary. The treatment arm = 1 when the patient is enrolled in the treatment group, treatment arm = 0 when in the control group.
}
  \item{prior.para}{
a list of hyper-parameters in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
}
}

\value{
a numeric vector; the posterior probabilities of the candidate models.
}
\references{
Duan, Y. and Parast, L., 2023. Flexible evaluation of surrogate markers with Bayesian model averaging. Statistics in Medicine.
}
\author{
Yunshan Duan
}
\examples{
data(exampleData)
post.model(Y = exampleData$Y, S = exampleData$S, A = exampleData$A)
}

