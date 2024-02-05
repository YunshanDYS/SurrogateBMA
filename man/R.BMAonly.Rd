\name{R.BMAonly}
\alias{R.BMAonly}
\title{
Calculates the proportion of treatment effect explained
}
\description{
Calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on the surrogate marker using Bayesian Model Averaging. This function is intended to be used for a fully observed continuous outcome and one single continuous surrogate marker. The user can also request to calculate a 95\% credible interval, evaluated by Bayesian bootstrapping. The MSE of predicting the primary outcome in the treatment can be calculated when test data is given.
}
\usage{
R.BMAonly(Y, S, A, nmc = 500, nBB = 100, conf.int = TRUE, alpha = 0.05, 
prior.para = NULL, testdata = FALSE, Ytest = NULL, Stest = NULL, Atest = NULL)
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
  \item{nmc}{
number of MCMC samples in posterior inference, default is 500.
}
  \item{nBB}{
number of replicates in Bayesian bootstrap, default is 100. Igored if conf.int = FALSE.
}
  \item{conf.int}{
TRUE or FALSE; indicates whether a 95\% credible interval for the proportion explained is requested, default is TRUE.
}
  \item{alpha}{
the confidence level for the credible interval, the 100(1 - \eqn{\alpha})\% credible interval is calculated. Default is 0.05, ignored if conf.int = FALSE.
}
  \item{prior.para}{
a list of hyper-parameters in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
}
\item{testdata}{
TRUE or FALSE; indicates whether test data is available, default is FALSE. The MSE of predicting the primary outcome is evaluated if TRUE. 
}
\item{Ytest}{
numeric vector; primary outcome in test data, assumed to be continuous. Default is NULL, ignored if testdata = FALSE.
}
  \item{Stest}{
numeric vector; surrogate marker in test data, assumed to be continuous. Default is NULL, ignored if testdata = FALSE.
}
  \item{Atest}{
numeric vector; treatment arm in test data, assumed to be binary. The treatment arm = 1 when the patient is enrolled in the treatment group, treatment arm = 0 when in the control group. Default is NULL, ignored if testdata = FALSE.
}
}

\value{
\item{R.BMA}{the estimate of the proportion explained by the surrogate marker.}
\item{p.model}{the posterior probability of the candidate models being true.}
\item{ci}{the credible interval for the proportion explained by the surrogate marker.}
\item{MSE}{the MSE of predicting the primary outcome for test data.}
}
\references{
Duan, Y. and Parast, L., 2023. Flexible evaluation of surrogate markers with Bayesian model averaging. Statistics in Medicine.
}
\author{
Yunshan Duan
}

\examples{
\donttest{
data(exampleData)
R.BMAonly(Y = exampleData$Y, S = exampleData$S, A = exampleData$A)
}
}


\keyword{internal}

