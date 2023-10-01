\name{R.BMAonly}
\alias{R.BMAonly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the proportion of treatment effect explained
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the proportion of treatment effect on the primary outcome explained by the treatment effect on the surrogate marker using Bayesian Model Averaging. This function is intended to be used for a fully observed continuous outcome and one single continuous surrogate marker. The user can also request to calculate a 95\% credible interval, evaluated by Bayesian bootstrapping. The MSE of predicting the primary outcome in the treatment can be calculated when test data is given.
}
\usage{
R.BMAonly(Y, S, A, nmc = 500, nBB = 100, conf.int = TRUE, alpha = 0.05, 
prior.para = NULL, testdata = FALSE, Ytest = NULL, Stest = NULL, Atest = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
numeric vector; primary outcome, assumed to be continuous.
}
  \item{S}{
%%     ~~Describe \code{S} here~~
numeric vector; surrogate marker, assumed to be continuous.
}
  \item{A}{
%%     ~~Describe \code{A} here~~
numeric vector; treatment arm, assumed to be binary. The treatment arm = 1 when the patient is enrolled in the treatment group, treatment arm = 0 when in the control group.
}
  \item{nmc}{
%%     ~~Describe \code{nmc} here~~
number of MCMC samples in posterior inference, default is 500.
}
  \item{nBB}{
%%     ~~Describe \code{nBB} here~~
number of replicates in Bayesian bootstrap, default is 100. Igored if conf.int = FALSE.
}
  \item{conf.int}{
%%     ~~Describe \code{conf.int} here~~
TRUE or FALSE; indicates whether a 95\% credible interval for the proportion explained is requested, default is TRUE.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
the confidence level for the credible interval, the 100(1 - \eqn{\alpha})\% credible interval is calculated. Default is 0.05, ignored if conf.int = FALSE.
}
  \item{prior.para}{
%%     ~~Describe \code{prior.para} here~~
a list of hyper-parameters in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
}
\item{testdata}{
%%     ~~Describe \code{testdata} here~~
TRUE or FALSE; indicates whether test data is available, default is FALSE. The MSE of predicting the primary outcome is evaluated if TRUE. 
}
\item{Ytest}{
%%     ~~Describe \code{Y} here~~
numeric vector; primary outcome in test data, assumed to be continuous. Default is NULL, ignored if testdata = FALSE.
}
  \item{Stest}{
%%     ~~Describe \code{S} here~~
numeric vector; surrogate marker in test data, assumed to be continuous. Default is NULL, ignored if testdata = FALSE.
}
  \item{Atest}{
%%     ~~Describe \code{A} here~~
numeric vector; treatment arm in test data, assumed to be binary. The treatment arm = 1 when the patient is enrolled in the treatment group, treatment arm = 0 when in the control group. Default is NULL, ignored if testdata = FALSE.
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{R.BMA}{the estimate of the proportion explained by the surrogate marker.}
\item{p.model}{the posterior probability of the candidate models being true.}
\item{ci}{the credible interval for the proportion explained by the surrogate marker.}
\item{MSE}{the MSE of predicting the primary outcome for test data.}
}
\references{
%% ~put references to the literature/web site here ~
Duan and Parast (2023+). Flexible Evaluation of Surrogate Markers with Bayesian Model Averaging.
}
\author{
%%  ~~who you are~~
Yunshan Duan
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
data(exampleData)
R.BMAonly(Y = exampleData$Y, S = exampleData$S, A = exampleData$A)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

\keyword{internal}

