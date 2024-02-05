\name{gen.prior}
\alias{gen.prior}
\title{
Generates the default prior hyperparameters
}
\description{
Generates a list of hyper-parameters we use by default in the inverse-Gamma-Normal prior for the variance and coefficients, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models. An Inv-Gamma(a0, b0) - Normal(mu0, \eqn{\sigma^2} Gamma0) prior is assumed.
}
\usage{
gen.prior()
}
\value{
A list of hyper-parameters, including a0_list, b0_list, mu0_list, Gamma0_list, Gamma0_inv_list , each being a list of 5 with 5 parameters under the 5 different candidate models.
}
\references{
Duan, Y. and Parast, L., 2023. Flexible evaluation of surrogate markers with Bayesian model averaging. Statistics in Medicine.
}
\author{
Yunshan Duan
}

\examples{
gen.prior()

# if we want to change the prior parameter a0 in model 1
prior.para <- gen.prior()
prior.para$a0_list[[1]] <- 1
}

