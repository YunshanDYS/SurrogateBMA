\name{post.theta}
\alias{post.theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generates posterior samples of the parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides a function to generate MCMC samples of p(theta_m | M = m, Data).
}
\usage{
post.theta(Y, S, A, m, nmc = 500, prior.para = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{nmc}{
%%     ~~Describe \code{nmc} here~~
}
  \item{prior.para}{
%%     ~~Describe \code{prior.para} here~~
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
a numeric matrix; each row contains a MC sample from the posterior of the parameters given a certain model.
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
%% \examples{
%% ##---- Should be DIRECTLY executable !! ----
%% ##-- ==>  Define data, use random,
%% ##--	or do  help(data=index)  for the standard data sets.

%% }
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
