# remotes::install_github("Ossifragus/OSMAC")
library(OSMAC)
library(MASS)

# this code is based on the algorithm for optimatl subsampling for logistic regression with balanced data
# described in "Optimal subsampling for large sample logistic regression" (Wang and others, 2018).
# Here we expand the original algorithm by choosing q based on power calculation for hypotheses testing.

# This function performs two-step optimal subsampling ang analysis, based on q0 (Step 1 sample size), required power,
# the coefficient that we want to test, and method (we use Wang and Others' original terms, "mmse" for A-optimal
# and "mvc" for L-optimal).

optimal_logistic_for_hypthoeses_testing <- function(X, Y, q0, power, alternative, coefficient, method) {
  n= nrow(X)
  n1 <- sum(Y)
  n0 <- n - n1
  PI.prop <- rep(1/(2 * n0), n)
  PI.prop[Y == 1] <- 1/(2 * n1)
  idx.prop <- sample(1:n, q0, T, PI.prop)
  x.prop <- X[idx.prop, ]
  y.prop <- Y[idx.prop]
  pinv.prop <- 1/PI.prop[idx.prop]
  fit.prop <- getMLE(x = x.prop, y = y.prop, w = pinv.prop)
  beta.prop <- fit.prop$par
  if (anyNA(beta.prop)) {
    result <- list(opt = NA, msg = "first stage not converge")
  }
  if (method == "mmse") {
    P.prop <- 1 - 1/(1 + exp(X %*% beta.prop))
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * 
                                     pinv.prop))
    PI.mMSE <- sqrt((Y - P.prop)^2 * rowSums((X %*% W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE) # These are the probabilities
    
    # estimate se based on Stage One:
    idx.prop <- sample(1:n, q0, T, PI.mMSE)
    x.prop <- X[idx.prop, ]
    y.prop <- Y[idx.prop]
    pinv.mMSE <- 1/PI.mMSE[idx.prop]
    p.mMSE <- 1 - 1/(1 + exp(c(x.prop %*% beta.prop)))
    w.mMSE <- p.mMSE * (1 - p.mMSE)
    W.mMSE <- solve(t(x.prop) %*% (x.prop * (w.mMSE * pinv.mMSE))) * q0 * n
    Vc.mMSE <- t(x.prop) %*% (x.prop * (y.prop - p.mMSE)^2 * pinv.mMSE^2)/q0^2/n^2
    V.mMSE <- W.mMSE %*% Vc.mMSE %*% W.mMSE * q0
    se <- sqrt(diag(V.mMSE))
    
    # calculate q
    q <- ceiling((qnorm(0.95)-qnorm(1-power))^2*se[coefficient]^2/(alternative^2))
    q <- min(n, q)
    
    # Stage Two
    idx.mMSE <- sample(1:n, q, T, PI.mMSE)
    x.mMSE <- X[c(idx.mMSE), ]
    y.mMSE <- Y[c(idx.mMSE)]
    pinv.mMSE <- 1/PI.mMSE[idx.mMSE]
    fit.mMSE <- getMLE(x = x.mMSE, y = y.mMSE, w = pinv.mMSE)
    beta.mMSE <- fit.mMSE$par
    p.mMSE <- 1 - 1/(1 + exp(c(x.mMSE %*% beta.mMSE)))
    w.mMSE <- p.mMSE * (1 - p.mMSE)
    W.mMSE <- solve(t(x.mMSE) %*% (x.mMSE * (w.mMSE * pinv.mMSE))) * q * n
    Vc.mMSE <- t(x.mMSE) %*% (x.mMSE * (y.mMSE - p.mMSE)^2 * pinv.mMSE^2)/q^2/n^2
    V.mMSE <- W.mMSE %*% Vc.mMSE %*% W.mMSE
    se <- sqrt(diag(V.mMSE))
    msg <- c(fit.prop$message, fit.mMSE$message)
    result <- list(par = beta.mMSE, se = se, msg = msg, 
                   method = "mmse", q=q, Z_score = beta.mMSE[coefficient]/se[coefficient])
  }
  
  else if(method == "mvc") {
    P.prop  <- 1 - 1 / (1 + exp(X %*% beta.prop))
    PI.mVc <- sqrt((Y - P.prop)^2 * rowSums(X^2))
    PI.mVc <- PI.mVc / sum(PI.mVc) # these are the weights
    
    # estimate se based on Stage 1
    idx.prop <- sample(1:n, q0, T, PI.mVc)
    x.prop <- X[idx.prop, ]
    y.prop <- Y[idx.prop]
    pinv.mVc <- 1/PI.mVc[idx.prop]
    p.mVc  <- 1 - 1 / (1 + exp(c(x.prop %*% beta.prop)))
    w.mVc <- p.mVc * (1 - p.mVc)
    W.mVc <- solve(t(x.prop) %*% (x.prop * (w.mVc * pinv.mVc))) * q0 * n
    Vc.mVc <- t(x.prop) %*% (x.prop * (y.prop-p.mVc)^2 * pinv.mVc^2) / q0^2 / n^2
    V.mVc <- W.mVc %*% Vc.mVc %*% W.mVc * q0
    se <- sqrt(diag(V.mVc))
    
    # calculate q
    q <- ceiling((qnorm(0.95)-qnorm(1-power))^2*se[coefficient]^2/(alternative^2))
    q <- min(q,n)
    
    # Stage Two
    idx.mVc <- sample(1:n, q, T, PI.mVc)
    x.mVc <- X[idx.mVc,]
    y.mVc <- Y[idx.mVc]
    pinv.mVc <- 1 / PI.mVc[idx.mVc]
    fit.mVc <- getMLE(x=x.mVc, y=y.mVc, w=pinv.mVc)
    
    beta.mVc <- fit.mVc$par
    p.mVc  <- 1 - 1 / (1 + exp(c(x.mVc %*% beta.mVc)))
    w.mVc <- p.mVc * (1 - p.mVc)
    W.mVc <- solve(t(x.mVc) %*% (x.mVc * (w.mVc * pinv.mVc))) * q * n
    Vc.mVc <- t(x.mVc) %*% (x.mVc * (y.mVc-p.mVc)^2 * pinv.mVc^2) / q^2 / n^2
    V.mVc <- W.mVc %*% Vc.mVc %*% W.mVc
    
    se <- sqrt(diag(V.mVc))
    msg <- c(fit.prop$message, fit.mVc$message)
    result <- list(par = beta.mVc, se = se, msg = msg, 
                   method = "mvc", q=q, Z_score = beta.mVc[coefficient]/se[coefficient])
  }
  
  result
}

# this function simulates data and perform analysis based on settings defined in Wang and Others (2018)
power_simulation <- function(n, q0, power, alternative, coefficient, r, method, distribution) {
  beta = rep(0.1, r)
  if(distribution == "mzNormal") {
    var.Matrix = matrix(0.5, nrow = r, ncol = r)
    diag(var.Matrix) = 1
    X = mvrnorm(n = n, mu = numeric(r), Sigma = var.Matrix)
    linear.Predictor = X %*% beta
    p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
    Y = rbinom(n = n, size = 1, prob = p)
  }
  else if(distribution == "nzNormal") {
    # nzNormal
    var.Matrix = matrix(0.5, nrow = r, ncol = r)
    diag(var.Matrix) = 1
    X = mvrnorm(n = n, mu = rep(1.5, r), Sigma = var.Matrix)
    linear.Predictor = X %*% beta
    p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
    Y = rbinom(n = n, size = 1, prob = p)
  }
  else if(distribution == "usNormal") {
    var.Matrix = matrix(0.5, nrow = r, ncol = r)
    diag(var.Matrix) = (1/(1:r))
    X = mvrnorm(n = n, mu = numeric(r), Sigma = var.Matrix, tol = 1)
    linear.Predictor = X %*% beta
    p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
    Y = rbinom(n = n, size = 1, prob = p)
  }
  else if(distribution == "T3") {
    var.Matrix = matrix(0.5, nrow = r, ncol = r)/10
    diag(var.Matrix) = 1/10
    X = rmvt(n = n, sigma = var.Matrix, df = 3)
    linear.Predictor = X %*% beta
    p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
    Y = rbinom(n = n, size = 1, prob = p)
  }
  else if(distribution == "exp") {
    X <- matrix(data = rexp(n = n * r, rate = 2), nrow = n, ncol = r)
    linear.Predictor = X %*% beta
    p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
    Y = rbinom(n = n, size = 1, prob = p)
  }
  
  optimal_logistic_for_hypthoeses_testing(X=X, Y=Y, q0, power, alternative, coefficient, method)
}

### run one simulations
power_simulation(n = 100000, q0 = 1000, power = 0.8, alternative = 0.1, coefficient = 7, r = 7, method = "mmse", distribution = "T3")

### verify that power is indead 80%:
mean(replicate(n=100,
               expr = power_simulation(n = 100000, q0 = 1000, power = 0.8, alternative = 0.1,
                                       coefficient = 7, r = 7, method = "mmse", distribution = "T3")$Z_score > qnorm(0.95)))