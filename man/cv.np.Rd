\name{cv.np}
\alias{cv.np}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate prodiction MSE for nonparametric method
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to calculate prodiction MSE for nonparametric method.
}
\usage{
cv.np(Y, S, A, Ytest, Stest, Atest)
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
  \item{Ytest}{
%%     ~~Describe \code{Y} here~~
}
  \item{Stest}{
%%     ~~Describe \code{S} here~~
}
  \item{Atest}{
%%     ~~Describe \code{A} here~~
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
a numeric value; MSE in cross-validation.
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
