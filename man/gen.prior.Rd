\name{gen.prior}
\alias{gen.prior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generates the default prior hyperparameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a list of hyper-parameters we use by default in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
}
\usage{
gen.prior()
}
%- maybe also 'usage' for other objects documented here.
%% \arguments{
%% }
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of hyper-parameters, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models.
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
gen.prior()

# if we want to change the prior parameter a0 in model 1
prior.para <- gen.prior()
prior.para$a0_list[[1]] <- 1
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
