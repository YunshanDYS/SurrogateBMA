\name{post.model}
\alias{post.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the posterior probability of the candidate models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Gives the posterior probability of each candidate model being true.
}
\usage{
post.model(Y, S, A, prior.para = NULL)
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
  \item{prior.para}{
%%     ~~Describe \code{prior.para} here~~
a list of hyper-parameters in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
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
a numeric vector; the posterior probabilities of the candidate models.
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
post.model(Y = exampleData$Y, S = exampleData$S, A = exampleData$A)
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
