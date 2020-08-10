
rm(list = ls())

library(BayesLogit)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)


source("simdata_biased_sample.R")
sourceCpp("sample_X.cpp")

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

# x, beta and gamma's
x_prec = diag(rep(0.001, d_x))
b_prec = diag(rep(0.001, 1 + d_x + d_u))
gam_prec = diag(rep(0.001, 1 + d_z + d_u))

PIF_rep = array(NA, dim = c(n_rep, 3, 5))  # true PIF, mean PIF, sd PIF, 2.5 quantile PIF, 97.5 quantile PIF

for (n_r in 1:n_rep) {
  
  # Generate data
  dataset = simdata(n = 10000)
  sum(dataset$data$R)  # phi_0=-2, n_v=966; phi_0=-2.8, n_v=483; phi_0=-3.5, n_v=273; phi_0=-4.5, n_v=105
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
  
  # storage objects
  BETA = matrix(NA, nr = n_length, nc = 1 + d_x + d_u)
  X_mis_str = array(NA, dim = c(n_length, n_m, d_x))
  GAMMA_str = array(NA, dim = c(n_length, d_x, 1 + d_z + d_u))
  INVSIGMA_str = array(NA, dim = c(n_length, d_x, d_x))
  PIF = array(NA, dim = c(n_length, 3))
  
  # initial values
  beta = as.numeric(coef(glm(Y ~ X.1 + X.2 + U, family = binomial, data=data[data$R == 1,])))
  X_obs = as.matrix(data[data$R == 1, c("X.1", "X.2")])
  mu_mis = colMeans(X_obs)
  sigma_mis = var(X_obs)
  X_mis = rep(1, n_m) %*% t(mu_mis) + 
    rmvnorm(n_m, mean = rep(0, d_x), sigma = sigma_mis)
  
  X_all_1 = cbind(1, as.matrix(rbind(X_mis, X_obs)))
  X_U_1 = cbind(X_all_1, U)
  
  Sigma = sigma_mis
  invSigma = solve(Sigma)
  
  gam1 = as.numeric(coef(lm(X.1 ~ Z.1 + Z.2 + U, data=data[data$R == 1,])))
  gam2 = as.numeric(coef(lm(X.2 ~ Z.1 + Z.2 + U, data=data[data$R == 1,])))
  Gamma = rbind(gam1, gam2)
  
  # wishart prior
  nu0 = d_x + 2
  S_0 = Sigma
  
  # Gibbs sampler
  for (i in 1:n_length) {
    
    # 1. sample latent variable
    eta = c(X_U_1 %*% beta)
    w = rpg(num = n, h = 1, z = eta)
    
    # 2. sample beta
    v_beta = solve(crossprod(X_U_1*w, X_U_1) + b_prec)
    m_beta = c(v_beta %*% (t(X_U_1) %*% (Y - 1/2)))
    beta = c(rmvnorm(1, m_beta, v_beta))
    BETA[i,] = beta
    
    # 3. sample unobserved x
    beta_nx = beta[-id_x]
    beta_x = beta[id_x]
    eta_nx = c(cbind(1, U) %*% beta_nx)
    zeta = Z_U_1 %*% t(Gamma)
    
    X_mis = sample_X(n_m, d_x, beta_x, tcrossprod(beta_x), w[1:n_m], eta_nx, zeta, invSigma, Y[1:n_m])
    X_all = as.matrix(rbind(X_mis, X_obs))
    X_all_1 = cbind(1, as.matrix(rbind(X_mis, X_obs)))
    X_U_1 = cbind(X_all_1, U)
    X_mis_str[i,,] = X_mis
    
    # 4. sample Gamma
    for(j in 1:nrow(Gamma)){
      v_gam_j = solve(invSigma[j,j] * crossprod(Z_U_1) + gam_prec)
      m_gam_j = v_gam_j %*% (crossprod(Z_U_1, X_all) %*% invSigma[,j] - 
                               crossprod(Z_U_1) %*% t(Gamma[-j,,drop=FALSE]) %*% invSigma[-j,j])
      Gamma[j,] = c(rmvnorm(1, m_gam_j, v_gam_j))
    }
    GAMMA_str[i,,] = Gamma
    
    # 5. sample Sigma
    SSR = crossprod(X_all - Z_U_1 %*% t(Gamma))
    invSigma = rwish(1, nu0 + n, solve(S_0 + SSR))
    INVSIGMA_str[i,,] = invSigma
    
    # 6. compute PIF
    X_U_1_star_x1 = X_U_1
    X_U_1_star_x2 = X_U_1
    X_U_1_star_x1x2 = X_U_1
    X_U_1_star_x1[,2] = X_U_1[,2] - 0.5
    X_U_1_star_x2[,3] = X_U_1[,3] - 0.5
    X_U_1_star_x1x2[,2] = X_U_1[,2] - 0.5
    X_U_1_star_x1x2[,3] = X_U_1[,3] - 0.5
    
    PIF[i,1] = 1 - mean(plogis(X_U_1_star_x1 %*% beta_0)) / mean(plogis(X_U_1 %*% beta_0))
    PIF[i,2] = 1 - mean(plogis(X_U_1_star_x2 %*% beta_0)) / mean(plogis(X_U_1 %*% beta_0))
    PIF[i,3] = 1 - mean(plogis(X_U_1_star_x1x2 %*% beta_0)) / mean(plogis(X_U_1 %*% beta_0))  
    
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
colMeans(GAMMA_str[-(1:n_burn),1,])
colMeans(GAMMA_str[-(1:n_burn),2,])

solve(colMeans(INVSIGMA_str[-(1:n_burn),,]))

traceplot(mcmc(PIF[,1]))
traceplot(mcmc(PIF[,2]))
traceplot(mcmc(PIF[,3]))
