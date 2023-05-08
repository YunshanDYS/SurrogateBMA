# Function to calculate F_m(theta_m)
R.m.theta <- function(S, m, theta, A, weight = FALSE, nBB = 100) {
  beta <- theta[2:5]
  sig2 <- theta[1]
  S1 <- S[which(A == 1)]
  S0 <- S[which(A == 0)]
  n1 <- length(S1)
  n0 <- length(S0)
  
  if (!weight) {
    if (m == 3){
      # Model 3: g(x) = x
      E1 <- mean((beta[1] + beta[2] + (beta[3] + beta[4])*S1) ) 
      E0 <- mean((beta[1] + beta[3]*S0) ) 
      AB <- mean((beta[1] + beta[2] + (beta[3] + beta[4])*S0) )
      Rvalue <- (E1 - AB)/(E1 - E0)
    } else if (m == 1) {
      # Model 1: g(x) = exp(x) - 1
      E1 <- mean((exp(beta[1] + beta[2] + (beta[3] + beta[4])*S1 + sig2/2) - 1) ) 
      E0 <- mean((exp(beta[1] + beta[3]*S0 + sig2/2) - 1) ) 
      AB <- mean((exp(beta[1] + beta[2] + (beta[3] + beta[4])*S0 + sig2/2) - 1) ) 
      Rvalue <- (E1 - AB)/(E1 - E0)
    } else if (m == 2) {
      # Model 2: g(x) = x^2
      E1 <- mean(((beta[1] + beta[2] + (beta[3] + beta[4])*S1)^2 + sig2) ) 
      E0 <- mean(((beta[1] + beta[3]*S0)^2 + sig2)) 
      AB <- mean(((beta[1] + beta[2] + (beta[3] + beta[4])*S0)^2 + sig2) ) 
      Rvalue <- (E1 - AB)/(E1 - E0)
    } else {
      # m = 4 or 5
      E1 <- gpxi_int(i = 1, m = m, beta = beta, S0 = S0, S1 = S1, sig2 = sig2)
      E0 <- gpxi_int(i = 0, m = m, beta = beta, S0 = S0, S1 = S1, sig2 = sig2)
      AB <- gpxi_int(i = 2, m = m, beta = beta, S0 = S0, S1 = S1, sig2 = sig2)
      Rvalue <- (E1 - AB)/(E1 - E0)
    }
  } else {
    if (m == 3){
      # Model 3: g(x) = x
      E1_vec <- (beta[1] + beta[2] + (beta[3] + beta[4]) * S1)
      E0_vec <- (beta[1] + beta[3] * S0)
      AB_vec <- (beta[1] + beta[2] + (beta[3] + beta[4]) * S0)
    } else if (m == 1) {
      # Model 1: g(x) = exp(x) - 1
      E1_vec <- (exp(beta[1] + beta[2] + (beta[3] + beta[4])*S1 + sig2/2) - 1)
      E0_vec <- (exp(beta[1] + beta[3]*S0 + sig2/2) - 1)
      AB_vec <- (exp(beta[1] + beta[2] + (beta[3] + beta[4])*S0 + sig2/2) - 1)
    } else if (m == 2) {
      # Model 2: g(x) = x^2
      E1_vec <- ((beta[1] + beta[2] + (beta[3] + beta[4])*S1)^2 + sig2) 
      E0_vec <- ((beta[1] + beta[3]*S0)^2 + sig2) 
      AB_vec <- ((beta[1] + beta[2] + (beta[3] + beta[4])*S0)^2 + sig2)
    } else {
      # m = 4 or 5
      E1_vec <- rep(NA, n1)
      for (k in 1:n1) {
        E1_vec[k] <- gpxi_int(i = 1, m = m, beta = beta, S0 = S0, S1 = c(S1[k]), sig2 = sig2)
      }
      E0_vec <- rep(NA, n0)
      for (k in 1:n0) {
        E0_vec[k] <- gpxi_int(i = 0, m = m, beta = beta, S0 = c(S0[k]), S1 = S1, sig2 = sig2)
      }
      AB_vec <- rep(NA, n0)
      for (k in 1:n0) {
        AB_vec[k] <- gpxi_int(i = 2, m = m, beta = beta, S0 = c(S0[k]), S1 = S1, sig2 = sig2)
      }
    }
    Rvalue <- c()
    for (b in 1:nBB) {
      # Bayesian bootstrap
      S1_samp <- sample(1:n1, n1, replace = T, prob = rdirichlet(n1, 1)) 
      S0_samp <- sample(1:n0, n0, replace = T, prob = rdirichlet(n0, 1)) 
      E1 <- mean(E1_vec[S1_samp]) 
      E0 <- mean(E0_vec[S0_samp]) 
      AB <- mean(AB_vec[S0_samp]) 
      R_temp <- (E1 - AB)/(E1 - E0)
      Rvalue <- c(Rvalue, R_temp)
    }
  }
  return(Rvalue)
}

gen.prior <- function() {
  M <- 5
  prior.para <- list()
  prior.para$a0_list <- prior.para$b0_list <- 
    prior.para$mu0_list <- prior.para$Gamma0_list <- 
    prior.para$Gamma0_inv_list <- list()
  for (m in 1:M) {
    prior.para$a0_list[[m]] <- 0.1
    prior.para$b0_list[[m]] <- 0.1
    prior.para$mu0_list[[m]] <- rep(0, 4)
    prior.para$Gamma0_list[[m]] <- diag(1 / 1000, 4)
    prior.para$Gamma0_inv_list[[m]] <- diag(1000, 4)
  }
  return(prior.para)
}

# Function to generate MC samples of p(theta_m \mid M = m, D)
post.theta <- function(Y, S, A, m, nmc = 500, prior.para = NULL) {
  M <- 5
  if (is.null(prior.para)) {
    # priors
    prior.para <- gen.prior()
  }
  a0 <- prior.para$a0_list[[m]]
  b0 <- prior.para$b0_list[[m]]
  mu0 <- prior.para$mu0_list[[m]]
  Gamma0 <- prior.para$Gamma0_list[[m]]
  n <- length(Y)
  # posterior for beta
  X <- cbind(rep(1,n), A, S, A*S)
  if (m == 1){
    # Model 1: g(x) = exp(x) - 1
    y <- log(Y + 1)
  }else if (m == 2){
    # Model 2: g(x) = x^2
    y <- sqrt(Y)
  }else if (m == 3){
    # Model 3: g(x) = x
    y <- Y
  }else if (m == 4){
    # Model 4: g(x) = \sqrt{x}
    y <- Y^2
  }else if (m == 5){
    # Model 5: g(x) = log(x + 1)
    y <- exp(Y) - 1
  }
  Gamman <- t(X)%*%X + Gamma0
  mun <- solve(Gamman) %*% (Gamma0%*%mu0 + t(X)%*%y) 
  an <- a0 + n/2
  bn <- b0 + 1/2*(sum(y^2) + t(mu0)%*%Gamma0%*%mu0 - t(mun)%*%Gamman%*%mun)
  sig2 <- rep(NA, nmc)
  beta <- matrix(NA, nrow = 4, ncol = nmc)
  for (i in 1:nmc){
    sig2[i] <- rinvgamma(1, an, bn)
    beta[, i] <- rmvnorm(1, mean = mun, sigma = sig2[i]*solve(Gamman))
  }
  theta.samples <- cbind(sig2, t(beta))
  return(theta.samples)
}


# Function to calculate p(M = m \mid D)
post.model <- function(Y, S, A, prior.para = NULL) {
  M <- 5 # num of candidate models
  if (is.null(prior.para)) {
    # priors
    prior.para <- gen.prior()
  }
  p.model <- rep(NA, M)
  for (m in 1:M){
    a0 <- prior.para$a0_list[[m]]
    b0 <- prior.para$b0_list[[m]]
    mu0 <- prior.para$mu0_list[[m]]
    Gamma0_inv <- prior.para$Gamma0_inv_list[[m]]
    n <- length(Y)
    X <- cbind(rep(1,n), A, S, A*S)
    if (m == 1){
      # Model 1: g(x) = exp(x) - 1
      y <- log(Y + 1)
      Jac <- 1/(Y + 1)
    }else if (m == 2){
      # Model 2: g(x) = x^2
      y <- sqrt(Y)
      Jac <- Y^(-0.5)/2
    }else if (m == 3){
      # Model 3: g(x) = x
      y <- Y
      Jac <- rep(1, n)
    }else if (m == 4){
      # Model 4: g(x) = \sqrt{x}
      y <- Y^2
      Jac <- 2*Y
    }else if (m == 5){
      # Model 5: g(x) = log(x + 1)
      y <- exp(Y) - 1
      Jac <- exp(Y)
    }
    if (is.infinite(sd(y))){
      p.model[m] <- -1e-10
    }else{
      # Marginal dist of Y \mid S, T, MVSt
      p.model[m] <- dmvt(y, delta = X%*%mu0, sigma = b0/a0*(diag(1,n) + X%*%Gamma0_inv%*%t(X)),
                       df = 2 * a0, log = TRUE) + sum(log(Jac))
    }
  }
  p.model[is.na(p.model)] <- -1e-10
  p.model <- p.model - max(p.model)
  p.model <- exp(p.model)
  p.model <- p.model/sum(p.model)
  return(p.model)
}

# Function R.BMA
R.BMA.est <- function(Y, S, A, 
                      nmc = 500, nBB = 100, conf.int = TRUE, alpha = 0.05,
                      prior.para =  NULL) {
  M <- 5 # num of candidate models
  # store output
  R_BMA_m <- rep(NA, M)
  p.model <- rep(NA, M)
  # calculate posterior prob of the models being true
  p.model <- post.model(Y, S, A, prior.para)
  # init
  n_samples_v <- round(p.model * nmc * nBB) # num of posterior samples in BB
  n <- length(A)
  n0 <- length(which(A == 0))
  n1 <- length(which(A == 1))
  Rsamples <- c()
  SA <- cbind(S, A)
  for (m in 1:M){
    R_BMA_post <- rep(NA, nmc)
    R_BMA_post_weighted <- c()
    if (p.model[m] < 1e-6){
      # avoid unnecessary calculation
      R_BMA_m[m] <- 0.5
    }else{
      theta_out <- post.theta(Y, S, A, m, nmc = nmc, prior.para = prior.para)
      for (i in 1:nmc){
        theta <- theta_out[i, ]
        R_BMA_post[i] <- R.m.theta(S, m, theta, A, weight = FALSE)
        if (conf.int) {
          R_temp <- R.m.theta(S, m, theta, A, weight = TRUE, nBB = nBB)
          R_BMA_post_weighted <- c(R_BMA_post_weighted, R_temp)
        }
      }
      R_BMA_m[m] <- mean(R_BMA_post)
      if (conf.int) {
        if (n_samples_v[m] != 0) {
          Rsamples <- c(Rsamples, R_BMA_post_weighted[sample(1:(nmc*nBB), n_samples_v[m])])
        }
      }
    }
  }
  R.BMA <- sum(p.model*R_BMA_m)
  ci <- NULL
  if (conf.int) {
    ci <- quantile(Rsamples, c((alpha / 2), (1 - alpha / 2)))
  }
  return(list(R.BMA = R.BMA, p.model = p.model, ci = ci))
}


