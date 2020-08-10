
rm(list = ls())

library(BayesLogit)
library(mvtnorm)

source("simdata.R")
# source("simdata_mvgamma.R")

# number of burn-ins and iterations
n_rep = 100
n_burn = 5000
n_iter = 5000
n_length = n_burn + n_iter

d_x = d_z = 2  # number of error-prone covariates
d_u = 1        # number of error-free covariates

id_x = 2:(1+d_x)

# true values for beta
beta_0 = c(-4, 1, 0.5, 0.5)  # decrease the default prevalence from -2 to -4
# beta_0 = c(-2, 1, 0.5, 0.5)

# prior parameters
b_prec = diag(rep(0.001, 1 + d_x + d_u))

PIF_rep = array(NA, dim = c(n_rep, 3, 5))   # true PIF, mean PIF, sd PIF, 2.5 quantile PIF, 97.5 quantile PIF

for (n_r in 1:n_rep) {
  
  # Generate data
  dataset = simdata(n_m = 10000, n_v = 100)
  data = dataset$data
  data_full = dataset$data_full
  n = nrow(data)
  n_v = sum(data$R)
  n_m = n - n_v
  U = data$U
  Z = cbind(data$Z.1, data$Z.2)
  Y = data$Y
  Z_U_1 = cbind(1, Z, U)
  
  # true PIF with half-unit decrease
  X_U_1_full = cbind(1, as.matrix(data_full[,c("X.1", "X.2", "U")]))
  X_U_1_full_star_x1x2 = X_U_1_full
  X_U_1_full_star_x1 = X_U_1_full
  X_U_1_full_star_x2 = X_U_1_full
  X_U_1_full_star_x1x2[,c("X.1", "X.2")] = X_U_1_full[,c("X.1", "X.2")] - 0.5
  X_U_1_full_star_x1[,"X.1"] = X_U_1_full[,"X.1"] - 0.5
  X_U_1_full_star_x2[,"X.2"] = X_U_1_full[,"X.2"] - 0.5
  
  PIF_rep[n_r,1,1] = 1 - mean(plogis(X_U_1_full_star_x1 %*% beta_0)) / mean(plogis(X_U_1_full %*% beta_0))
  PIF_rep[n_r,2,1] = 1 - mean(plogis(X_U_1_full_star_x2 %*% beta_0)) / mean(plogis(X_U_1_full %*% beta_0))
  PIF_rep[n_r,3,1] = 1 - mean(plogis(X_U_1_full_star_x1x2 %*% beta_0)) / mean(plogis(X_U_1_full %*% beta_0))
  
  Z_U_1_star_x1x2 = Z_U_1
  Z_U_1_star_x1x2[,id_x] = Z_U_1[,id_x] - 0.5
  Z_U_1_star_x1 = Z_U_1
  Z_U_1_star_x1[,2] = Z_U_1[,2] - 0.5
  Z_U_1_star_x2 = Z_U_1
  Z_U_1_star_x2[,3] = Z_U_1[,3] - 0.5
  
  # storage objects
  BETA = matrix(NA, nr = n_length, nc = 1 + d_x + d_u)
  PIF = array(NA, dim = c(n_length, 3))
  
  beta = as.numeric(coef(glm(Y ~ Z.1 + Z.2 + U, family = binomial, data=data)))
  
  # Gibbs sampler
  for (i in 1:n_length) {
    
    # 1. sample latent variable
    eta = c(Z_U_1 %*% beta)
    w = rpg(num = n, h = 1, z = eta)
    
    # 2. sample beta
    v_beta = solve(crossprod(Z_U_1*w, Z_U_1) + b_prec)
    m_beta = c(v_beta %*% (t(Z_U_1) %*% (Y - 1/2)))
    beta = c(rmvnorm(1, m_beta, v_beta))
    BETA[i,] = beta
    
    # 3. compute PIF
    PIF[i,1] = 1 - mean(plogis(Z_U_1_star_x1x2 %*% beta)) / mean(plogis(Z_U_1 %*% beta))
    PIF[i,2] = 1 - mean(plogis(Z_U_1_star_x1 %*% beta)) / mean(plogis(Z_U_1 %*% beta))
    PIF[i,3] = 1 - mean(plogis(Z_U_1_star_x2 %*% beta)) / mean(plogis(Z_U_1 %*% beta))
    
    cat("\r", paste("repetition", n_r, "iteration", i, sep=" "))
    flush.console()
  }
  
  PIF_rep[n_r,,2] = apply(PIF[-(1:n_burn),], 2, mean)
  PIF_rep[n_r,,3] = apply(PIF[-(1:n_burn),], 2, sd)
  PIF_rep[n_r,,4] = apply(PIF[-(1:n_burn),], 2, function(x) quantile(x, 0.025))
  PIF_rep[n_r,,5] = apply(PIF[-(1:n_burn),], 2, function(x) quantile(x, 0.975))
}

print('----------------------------------------------------------------------')

print('PIF bias percentage')
apply(PIF_rep[,,2] - PIF_rep[,,1], 2, mean) / apply(PIF_rep[,,1], 2, mean)

print('PIF coverage')
apply((PIF_rep[,,4] <= PIF_rep[,,1])*(PIF_rep[,,5] >= PIF_rep[,,1]), 2, mean)

print('PIF MSE')
apply((PIF_rep[,,2] - PIF_rep[,,1])^2, 2, mean)

print('PIF ASE')
apply(PIF_rep[,,3], 2, mean)

print('PIF MCSE')
apply(PIF_rep[,,2], 2, sd)


library(coda)
traceplot(mcmc(BETA[,1]))
traceplot(mcmc(BETA[,2]))
traceplot(mcmc(BETA[,3]))
traceplot(mcmc(BETA[,4]))

colMeans(BETA[-(1:n_burn),])

traceplot(mcmc(PIF[,1]))
traceplot(mcmc(PIF[,2]))
traceplot(mcmc(PIF[,3]))
