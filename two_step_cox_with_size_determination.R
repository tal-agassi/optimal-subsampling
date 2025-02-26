library(survival)
library(MASS)
library(Rcpp)
library(Epi)
library(multipleNCC)
library(msm)
library(dplyr)
sourceCpp("cox-subsampling-expansion.cpp", verbose=TRUE)
library(tidyverse)
library(ggplot2)

information_score_matrix <- function(beta,weights=NULL,times,truncs=NULL,status,covariates,samp_prob=NULL,information_mat=T,score_res_mat=T,diagon=F) {
  if(is.null(weights)) {weights = c(-1)}
  if(is.null(truncs)) {truncs = c(-1)}
  if(is.null(samp_prob))
  {
    ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1*information_mat,1*score_res_mat,1*diagon)
    if(information_mat)
    {
      tmp = ret$hess  
      ret_hess = tmp + t(tmp)
      diag(ret_hess) = diag(tmp)
      ret[["hess"]] = -ret_hess 
    }
    return(ret)
  }
  if(samp_prob == "A")
  {
    ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1,1,0)
    tmp = ret$hess  
    ret_hess = tmp + t(tmp)
    diag(ret_hess) = diag(tmp)
    ret[["hess"]] = -ret_hess 
    ret[["samp_prob"]] = rcpp_A_OPT(ret[["residual"]],solve(ret[["hess"]]))
    return(ret)
  }
  if(samp_prob == "L")
  {
    if(information_mat)
    {
      ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,1,1,0)
      tmp = ret$hess  
      ret_hess = tmp + t(tmp)
      diag(ret_hess) = diag(tmp)
      ret[["hess"]] = -ret_hess 
    } else
    {
      ret = rcpp_score_wrapper(beta,weights,times,truncs,status,covariates,0,1,0)
    }
    ret[["samp_prob"]] = rcpp_L_OPT(ret[["residual"]])
    return(ret)
  }
}

two_step_cox_with_size_determination <- function(v, delta, X, q0_controls, alternative, power, coefficient, method) {
  
  n_events = sum(delta)
  n_cens <- n - n_events
  q0 = round(n_events*q0_controls)
  
  cens_ind = which(delta == 0)
  
  
  unif_cont_ind = sample(cens_ind,q0,replace = T) # uniform sampling
  samp_ind_unif = c(unif_cont_ind,setdiff(1:n,cens_ind))
  cens_weights_unif = length(cens_ind)/ q0
  weights_unif = ifelse(delta,1,cens_weights_unif)
  fit_samp_unif = coxph(Surv(time=v[samp_ind_unif],event=delta[samp_ind_unif]) ~ X[samp_ind_unif,],weights = weights_unif[samp_ind_unif],robust = F)  #no truncation
  unif_res <- coef(fit_samp_unif)
  
  # A optimality weights calculation + estimation of information matrix
  if(method=='A') {
    tmp_A1 = information_score_matrix(unif_res[m,],times = v,status = delta,covariates = X, samp_prob = "A",information_mat = F) #no truncation
    delta_A = delta[tmp_A1$ord]
    cens_ind_ord = which(!delta_A)
    
    # random sampling with score-optimal weights among censored
    samp_ind_cens = sample(1:n_cens,q0,replace = T,prob = tmp_A1$samp_prob)
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_A))
    cens_weights_opt = (1/(tmp_A1$samp_prob * q0))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_A1$ord[samp_ind_opt]
    
    fit_samp_opt_A = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
    opt_A_res <- coef(fit_samp_opt_A)
    
    tmp_A2 = information_score_matrix(opt_A_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
    ind_rm = (1:n_events)+q0
    order_rm = tmp_A2$ord[!(tmp_A2$ord %in% ind_rm)] 
    Score_A = tmp_A2$residual / tmp_A1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_A = cov(Score_A)
    I_inv_A = solve(tmp_A2$hess)
    var_opt_A = I_inv_A + I_inv_A %*% phi_mat_A %*% I_inv_A/q0
    
    # calculating q
    
    q <- ceiling(((I_inv_A %*% phi_mat_A %*% I_inv_A)[coefficient,coefficient]*(qnorm(0.95)-qnorm(1-power))^2)/(alternative^2-I_inv_A[coefficient,coefficient]*(qnorm(0.95)-qnorm(1-power))^2))
    
    # re-calculating A with the new q:
    
    # random sampling with score-optimal weights among censored
    samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_A1$samp_prob)
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_A))
    cens_weights_opt = (1/(tmp_A1$samp_prob * q))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_A1$ord[samp_ind_opt]
    
    fit_samp_opt_A = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
    opt_A_res <- coef(fit_samp_opt_A)
    
    tmp_A2 = information_score_matrix(opt_A_res,weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
    ind_rm = (1:n_events)+q
    order_rm = tmp_A2$ord[!(tmp_A2$ord %in% ind_rm)] 
    Score_A = tmp_A2$residual / tmp_A1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_A = cov(Score_A)
    I_inv_A = solve(tmp_A2$hess)
    var_opt_A = I_inv_A + I_inv_A %*% phi_mat_A %*% I_inv_A/q
    
    Z_score <- opt_A_res[coefficient]/sqrt(var_opt_A[coefficient,coefficient])
    
    results = list("coef" = opt_A_res, "var" = var_opt_A, "q"=q)
    
  }
  
  if(method=='L') {
    # L optimality weights calculation + estimation of information matrix
    
    tmp_L1 = information_score_matrix(unif_res,times = v,status = delta,covariates = X, samp_prob = "L",information_mat = F) #no truncation
    delta_L = delta[tmp_L1$ord]
    cens_ind_ord = which(!delta_L)
    
    # random sampling with score-optimal weights among censored
    samp_ind_cens = sample(1:n_cens,q0,replace = T,prob = tmp_L1$samp_prob)
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_L))
    cens_weights_opt = (1/(tmp_L1$samp_prob * q0))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_L1$ord[samp_ind_opt]
    
    fit_samp_opt_L = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
    opt_L_res <- coef(fit_samp_opt_L)
    
    tmp_L2 = information_score_matrix(opt_L_res,weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
    ind_rm = (1:n_events)+q0
    order_rm = tmp_L2$ord[!(tmp_L2$ord %in% ind_rm)]
    Score_L = tmp_L2$residual / tmp_L1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_L = cov(Score_L)
    I_inv_L = solve(tmp_L2$hess)
    var_opt_L = I_inv_L + I_inv_L %*% phi_mat_L %*% I_inv_L/q0
    
    # calculating q
    
    q <- ceiling(((I_inv_L %*% phi_mat_L %*% I_inv_L)[coefficient,coefficient]*(qnorm(0.95)-qnorm(1-power))^2)/(alternative^2-I_inv_L[coefficient,coefficient]*(qnorm(0.95)-qnorm(1-power))^2))
    
    # re-calculating L with the new q:
    
    # random sampling with score-optimal weights among censored
    samp_ind_cens = sample(1:n_cens,q,replace = T,prob = tmp_L1$samp_prob)
    samp_ind_opt = c(cens_ind_ord[samp_ind_cens],which(delta_L))
    cens_weights_opt = (1/(tmp_L1$samp_prob * q))[samp_ind_cens]
    weights_opt = c(cens_weights_opt,rep(1,n_events))
    
    samp_ord = tmp_L1$ord[samp_ind_opt]
    
    fit_samp_opt_L = coxph(Surv(time=v[samp_ord],event=delta[samp_ord]) ~ X[samp_ord,],weights = weights_opt,robust = F,init = unif_res[m,])  #no truncation
    opt_L_res <- coef(fit_samp_opt_L)
    
    tmp_L2 = information_score_matrix(opt_L_res[m,],weights = weights_opt,times = v[samp_ord],status = delta[samp_ord],covariates = X[samp_ord,]) #no truncation
    ind_rm = (1:n_events)+q
    order_rm = tmp_L2$ord[!(tmp_L2$ord %in% ind_rm)]
    Score_L = tmp_L2$residual / tmp_L1$samp_prob[samp_ind_cens][order_rm]
    phi_mat_L = cov(Score_L)
    I_inv_L = solve(tmp_L2$hess)
    var_opt_L = I_inv_L + I_inv_L %*% phi_mat_L %*% I_inv_L/q
    
    results = list("coef" = opt_L_res, "var" = var_opt_L, "q"=q)
  }
  
}


### simulating data
n = 100000 # sample size
r = 6 # number of coefficients
b <- c(3,-5,1,-1,1,-3,rep_len(c(-1,1/2),r-6))/10
upper_unif =  c(4,4,4,4,4,4,rep_len(c(1,1),r-6))
rates = c(0.001,0.005) 
c = rexp(n,0.2) # censoring times
X_vec = runif(n*r,0,rep(upper_unif,each = n))
X = matrix(X_vec,nrow = n,ncol = r)

linear_comb = X %*% b

knots = c(0,6)
obs_rates = exp(linear_comb) %*% rates
y = apply(obs_rates,1,rpexp,n=1,t=knots) 
v = pmin(c,y) # observed time
delta = v == y

write_csv(data.frame(v=v, delta=delta, X=X), 'C:/Users/טל/OneDrive - mail.tau.ac.il/שולחן העבודה/thesis/code for git/two_step_cox_data_example.csv')
### run on simulated data

two_step_cox_with_size_determination(v = v, delta = delta, X = X, q0_controls = 1, alternative = 0.3, power = 0.8, coefficient = 4, method = 'A')