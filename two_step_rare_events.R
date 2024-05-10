# remotes::install_github("Ossifragus/OSMAC")
library(OSMAC)
library(MASS)

### This function performs the two-step algorithm for logistic regression with rare events, which
### keeps all the events in the subsample. It doesn't include step 1.5.


twostep_rare_events <- function(X, Y, q0, q, method) {
  n0_ind <- which(Y==0)
  n1_ind <- which(Y==1)
  n_events <- length(n1_ind)
  unif_ind <- c(n1_ind, sample(n0_ind, q0, T))
  X_unif <- X[unif_ind,]
  Y_unif <- Y[unif_ind]
  cens_weights_unif <- length(n0_ind) / q0
  weights_unif = ifelse(Y[unif_ind]==1,1,cens_weights_unif)
  fit_unif <- getMLE(x = X_unif, y = Y_unif, w = weights_unif)
  beta_unif <- fit_unif$par
  
  if(method=='A') {
    
    P.unif  <- 1 - 1 / (1 + exp(X %*% beta_unif))
    p.unif <- P.unif[unif_ind]
    w.unif <- p.unif * (1 - p.unif)
    W.unif <- solve(t(X_unif) %*% (X_unif * w.unif * weights_unif))
    PI.mMSE <- sqrt((Y - P.unif)^2 * rowSums((X%*%W.unif)^2))
    PI.mMSE[n0_ind] <- PI.mMSE[n0_ind] / sum(PI.mMSE[n0_ind])
    idx.mMSE <- sample(n0_ind, q, T, PI.mMSE[n0_ind])
    x.mMSE <- X[c(idx.mMSE, n1_ind),]
    y.mMSE <- Y[c(idx.mMSE, n1_ind)]
    pinv.mMSE <- c(1 / (PI.mMSE[idx.mMSE] * q), rep(1, length(n1_ind)))
    fit.mMSE <- glm(y.mMSE~x.mMSE[,-1], family = "binomial", weights = pinv.mMSE)
    
    p.opt <- as.vector((1 - 1 / (1 + exp(X[c(idx.mMSE, n1_ind),] %*% fit.mMSE$coefficients))))
    w.opt <- p.opt * (1 - p.opt) 
    inv.M <- solve(t(x.mMSE) %*% (x.mMSE * (w.opt * pinv.mMSE)))
    
    p.opt <- as.vector((1 - 1 / (1 + exp(X[idx.mMSE,] %*% fit.mMSE$coefficients))))
    phi <- cov(x.mMSE[1:q,] * p.opt * pinv.mMSE[1:q] * q)
    
    var.est <- inv.M + inv.M %*% phi %*% inv.M / q
    se <- sqrt(diag(var.est))
    
    result <- list(par = fit.mMSE$coefficients, se = se, 
                   method = "mmse", q0=q0, q=q)
  }
  
  
  if(method == 'L') {
    
    P.unif <- 1 - 1/(1 + exp(X %*% beta_unif))
    PI.mVc <- sqrt((Y - P.unif)^2 * rowSums(X^2))
    PI.mVc[n0_ind] <- PI.mVc[n0_ind]/sum(PI.mVc[n0_ind])
    idx.mVc <- sample(n0_ind, q, T, PI.mVc[n0_ind])
    x.mVc <- X[c(idx.mVc, n1_ind), ]
    y.mVc <- Y[c(idx.mVc, n1_ind)]
    pinv.mVc <- c(1/(PI.mVc[idx.mVc] * q), rep(1, length(n1_ind)))
    fit.mVc <- glm(y.mVc~x.mVc[,-1], family = "binomial", weights = pinv.mVc)
    
    p.opt <- (1 - 1 / (1 + exp(X %*% fit.mVc$coefficients)))[c(idx.mVc, n1_ind)]
    w.opt <- p.opt * (1 - p.opt) 
    inv.M <- solve(t(x.mVc) %*% (x.mVc * (w.opt * pinv.mVc)))
    
    p.opt <- p.opt[1:length(idx.mVc)]
    phi <- cov(x.mVc[1:q,] * p.opt * pinv.mVc[1:q] * q)
    
    var.est <- inv.M + inv.M %*% phi %*% inv.M / q
    se <- sqrt(diag(var.est))
    
    result <- list(par =  fit.mVc$coefficients, se = se, 
                   method = "L", q0=q0, q=q)
  }
  
  result

}

### simulate example data
n=10000
r = 6
beta = rep(0.1, r)
var.Matrix = matrix(0.5, nrow = r, ncol = r)
diag(var.Matrix) = 1
X = mvrnorm(n = n, mu = numeric(r), Sigma = var.Matrix)
linear.Predictor = X %*% beta
p = exp(linear.Predictor)/(exp(linear.Predictor)+1)
Y = rbinom(n = n, size = 1, prob = p)

q0 = 1000
q = 1000
twostep_rare_events(X = X, Y = Y, q0 = q0, q = q, method = 'A')
twostep_rare_events(X = X, Y = Y, q0 = q0, q = q, method = 'L')
