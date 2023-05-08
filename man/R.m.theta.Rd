\name{R.m.theta}
\alias{R.m.theta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the R value given model and parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides a function to calculate F_m(theta_m).
}
\usage{
R.m.theta(S, m, theta, A, weight = FALSE, nBB = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{nBB}{
%%     ~~Describe \code{nBB} here~~
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
a numeric value or a numeric verctor when weight = TRUE. the value of R or the values of R using Bayesian bootstap replicates of the surrogate marker.
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
