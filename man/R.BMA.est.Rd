\name{R.BMA.est}
\alias{R.BMA.est}
\title{
Calculates the proportion of treatment effect explained
}
\description{
Calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on the surrogate marker using Bayesian Model Averaging (BMA). This function is intended to be used for a fully observed continuous outcome and one single continuous surrogate marker. The user can also request to calculate a 95\% credible interval, evaluated by Bayesian bootstrapping. A robust option is provided where we choose a more preferrable approach between the BMA and non-parametric methods based on cross-validation. 
}
\usage{
R.BMA.est(Y, S, A, method = "BMA", nmc = 500, nBB = 100, conf.int = TRUE, 
alpha = 0.05, prior.para = NULL, kfold.k = 3)
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
  \item{method}{
"BMA" or "robust"; Default is "BMA". If method = "robust", we choose a more preferrable approach between the BMA and non-parametric methods based on cross-validation. 
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
  \item{kfold.k}{
the value of k in k-fold cross validation, default is 3. Ignored if method = "BMA".
}
}

\value{
\item{R.est}{the estimate of the proportion explained by the surrogate marker.}
\item{p.model}{the posterior probability of the candidate models being true.}
\item{ci}{the credible (confidence) interval for the proportion explained by the surrogate marker.}
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
R.BMA.est(Y = exampleData$Y, S = exampleData$S, A = exampleData$A)
}
}

