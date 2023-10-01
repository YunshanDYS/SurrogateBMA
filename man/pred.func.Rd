\name{pred.func}
\alias{pred.func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate the expected primary outcome in the treatment group given the model and the parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides a function to calculate the expected primary outcome in the treatment group given the model and the parameters.
}
\usage{
pred.func(theta, ST, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{Y} here~~
}
  \item{ST}{
%%     ~~Describe \code{S} here~~
}
  \item{m}{
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
a numeric vector.
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
