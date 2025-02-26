# remotes::install_github("Ossifragus/OSMAC")
library(OSMAC)
library(MASS)
library(mvtnorm)
library(tidyverse)


two_step_rare_events_with_size_determination = function(x,y,method, coefficient, alternative, power) {
  n0_ind <- which(Y==0)
  n1_ind <- which(Y==1)
  n1 <- length(n1_ind)
  q0 = length(n1_ind)
  unif_ind <- c(n1_ind, sample(n0_ind, q0, T))
  X_unif <- X[unif_ind,]
  Y_unif <- Y[unif_ind]
  cens_weights_unif <- length(n0_ind)/q0
  weights_unif = ifelse(Y==1,1,cens_weights_unif)[unif_ind]
  fit_unif <- getMLE(x = X_unif, y = Y_unif, w = weights_unif)
  beta_unif <- fit_unif$par
  
  P.prop.S1  <- 1 - 1 / (1 + exp(X %*% beta_unif))
  p.prop.S1 <- P.prop.S1[unif_ind]
  w.prop.S1 <- p.prop.S1 * (1 - p.prop.S1)
  W.prop.S1 <- solve(t(X_unif) %*% (X_unif * w.prop.S1 * weights_unif))
  
  ### A opt
  if(method=='A') {
    # estimating se and phi for calculation of q

    PI.mMSE.S1 <- sqrt((Y - P.prop.S1)^2 * rowSums((X%*%W.prop.S1)^2))
    PI.mMSE.S1[n0_ind] <- PI.mMSE.S1[n0_ind] / sum(PI.mMSE.S1[n0_ind])
    idx.mMSE.S1 <- sample(n0_ind, q0, T, PI.mMSE.S1[n0_ind])
    x.mMSE.S1 <- X[c(idx.mMSE.S1, n1_ind),]
    pinv.mMSE.S1 <- c(1 / (PI.mMSE.S1[idx.mMSE.S1] * q0), rep(1, length(n1_ind)))
    p.opt.S1 <- as.vector((1 - 1 / (1 + exp(X[c(idx.mMSE.S1, n1_ind),] %*% beta_unif))))
    w.opt.S1 <- p.opt.S1 * (1 - p.opt.S1) 
    inv.M <- solve(t(x.mMSE.S1) %*% (x.mMSE.S1 * (w.opt.S1 * pinv.mMSE.S1)))
    p.opt.S1 <- p.opt.S1[1:length(idx.mMSE.S1)]
    phi <- cov(x.mMSE.S1[1:q0,] * p.opt.S1 * pinv.mMSE.S1[1:q0] * q0)
    
    q <- ceiling((inv.M %*% phi %*% inv.M)[coefficient,coefficient] * (qnorm(0.95) - qnorm(1-power))^2 / (alternative^2 - inv.M[coefficient,coefficient] * (qnorm(0.95)-qnorm(1-power))^2))
    
    # Stage 2
    idx.mMSE <- sample(n0_ind, q, T, PI.mMSE.S1[n0_ind])
    x.mMSE <- X[c(idx.mMSE, n1_ind),]
    y.mMSE <- Y[c(idx.mMSE, n1_ind)]
    pinv.mMSE <- c(1 / (PI.mMSE.S1[idx.mMSE] * q), rep(1, length(n1_ind)))
    fit.mMSE <- glm(y.mMSE~x.mMSE[,-1], family = "binomial", weights = pinv.mMSE)
    
    p.opt <- as.vector((1 - 1 / (1 + exp(X[c(idx.mMSE, n1_ind),] %*% fit.mMSE$coefficients))))
    w.opt <- p.opt * (1 - p.opt) 
    inv.M <- solve(t(x.mMSE) %*% (x.mMSE * (w.opt * pinv.mMSE)))
    
    p.opt <- p.opt[1:length(idx.mMSE)]
    phi <- cov(x.mMSE[1:q,] * p.opt * pinv.mMSE[1:q] * q)
    
    var.est <- inv.M + inv.M %*% phi %*% inv.M / q
    
    results = list("coefficients" = fit.mMSE$coefficients, 'sd' = sqrt(diag(var.est)), 'q'=q, 'q0'=q0, 'method'=method)
  }
  
  
  ### L opt
  
  if(method=='L')
    
  {  # estimating se and phi for calculation of q
    PI.mVc.S1 <- sqrt((Y - P.prop.S1)^2 * rowSums(X^2))
    PI.mVc.S1[n0_ind] <- PI.mVc.S1[n0_ind]/sum(PI.mVc.S1[n0_ind])
    idx.mVc.S1 <- sample(n0_ind, q0, T, PI.mVc.S1[n0_ind])
    x.mVc.S1 <- X[c(idx.mVc.S1, n1_ind), ]
    pinv.mVc.S1 <- c(1/(PI.mVc.S1[idx.mVc.S1] * q0), rep(1, length(n1_ind)))
    p.opt.S1 <- as.vector((1 - 1 / (1 + exp(X[c(idx.mVc.S1, n1_ind),] %*% beta_unif))))
    w.opt.S1 <- p.opt.S1 * (1 - p.opt.S1)
    inv.M <- solve(t(x.mVc.S1) %*% (x.mVc.S1 * (w.opt.S1 * pinv.mVc.S1)))
    
    p.opt.S1 <- p.opt.S1[1:length(idx.mVc.S1)]
    phi <- cov(x.mVc.S1[1:q0,] * p.opt.S1 * pinv.mVc.S1[1:q0] * q0)
    
    var.est <- inv.M + inv.M %*% phi %*% inv.M / q0
    
    q <- ceiling((inv.M %*% phi %*% inv.M)[coefficient,coefficient] * (qnorm(0.95) - qnorm(1-power))^2 / (alternative^2 - inv.M[coefficient,coefficient] * (qnorm(0.95)-qnorm(1-power))^2))
    
    #Stage Two
    
    idx.mVc <- sample(n0_ind, q, T, PI.mVc.S1[n0_ind])
    x.mVc <- X[c(idx.mVc, n1_ind), ]
    y.mVc <- Y[c(idx.mVc, n1_ind)]
    pinv.mVc <- c(1/(PI.mVc.S1[idx.mVc] * q), rep(1, length(n1_ind)))
    fit.mVc <- glm(y.mVc~x.mVc[,-1], family = "binomial", weights = pinv.mVc)
    
    p.opt <- as.vector((1 - 1 / (1 + exp(X[c(idx.mVc, n1_ind),] %*% fit.mVc$coefficients))))
    w.opt <- p.opt * (1 - p.opt) 
    inv.M <- solve(t(x.mVc) %*% (x.mVc * (w.opt * pinv.mVc)))
    
    p.opt <- p.opt[1:length(idx.mVc)]
    phi <- cov(x.mVc[1:q,] * p.opt * pinv.mVc[1:q] * q)
    
    var.est <- inv.M + inv.M %*% phi %*% inv.M / q
   
    results = list("coefficients" = fit.mVc$coefficients, 'sd' = sqrt(diag(var.est)), 'q'=q, 'q0'=q0, method'=method) 
  }
  
  results

}


### simulate data

r = 6 # num of coefficients
N = 100000 # sample size
beta = c(-4, rep(0.1, r)) # include intercept
var.Matrix = matrix(0.5, nrow = r, ncol = r)
diag(var.Matrix) = 1
X = cbind(rep(1, N), mvrnorm(n = N, mu = numeric(r), Sigma = var.Matrix))
linear.Predictor = X %*% beta
p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
Y = rbinom(n = N, size = 1, prob = p)

### run on simulated data
two_step_rare_events_with_size_determination(x = X, y = Y,method = 'A', coefficient = 5, alternative = 0.15, power = 0.8)
two_step_rare_events_with_size_determination(x = X, y = Y,method = 'L', coefficient = 5, alternative = 0.1, power = 0.85)
