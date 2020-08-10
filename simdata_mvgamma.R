### data generation
### Start Simple!
### basic scenario: two x's and one u

require(mvtnorm)

rmvgamma = function(n, shape=1, rate=1, corr=diag(length(shape))) {
  ## extract parameters, do sanity checks, deal with univariate case
  
  if(!is.matrix(corr) || !isSymmetric(corr))
    stop("'corr' must be a symmetric matrix")
  D = ncol(corr)
  
  Ds = length(shape)
  if(Ds > D)
    warning("'shape' longer than width of 'corr', truncating to fit")
  if(Ds != D)
    shape = rep(shape, length.out=D)
  
  Dr = length(rate)
  if(Dr > D)
    warning("'rate' longer than width of 'corr', truncating to fit")
  if(Dr != D)
    rate = rep(rate, length.out=D)
  
  if(D == 1) rgamma(n, shape, rate)
  
  ## generate standard multivariate normal matrix, convert to CDF
  
  Z = rmvnorm(n, sigma=corr)
  cdf = pnorm(Z)
  
  ## convert to gamma, return
  
  sapply(1:D, function(d) qgamma(cdf[,d], shape[d], rate[d]))
}

simdata = function(n_m, n_v){
  
  n = n_m + n_v  # main study and validation study sizes
  d_x = d_z = 2  # number of error-prone covariates
  d_u = 1        # number of error-free covariates
  
  # true values for beta
  beta_0 = c(-4, 1, 0.5, 0.5) 
  
  # true values for alpha, Gamma, Sigma_x
  Gamma_0 = matrix(c(-0.2, 0.6, 0.2, -0.7,  # moderate measurement error (correlation between X and Z)
                     0.1, 0.2, 0.4, -0.5), nr = d_x, nc = 1+d_x+d_u, byrow = TRUE)
  Sigma_x0 = matrix(c(1, 0.2,
                      0.2, 1), nr = d_x, nc = d_x)
  
  # true values for mu_z, Sigma_z
  mu_z0 = c(-0.2, 0.8) 
  Sigma_z0 = matrix(c(1, 0.2,
                      0.2, 1), nr = d_z, nc = d_z)
  
  # generate U and Z
  U = matrix(rnorm(n*d_u), nr = n, nc = d_u)
  Z = rmvnorm(n = n, mean = mu_z0, sigma = Sigma_z0)
  
  # generate X
  mu_x0 = cbind(1, Z, U) %*% t(Gamma_0)
  X = mu_x0 + rmvgamma(n = n, shape = 2, rate = 2, corr = Sigma_x0) - c(1, 1)
  
  # generate Y
  prob_Y = plogis(cbind(1, X, U) %*% beta_0)
  Y = rbinom(n, 1, prob_Y)
  
  # observed data
  R = rep(c(0, 1), c(n_m, n_v))
  data_full = data.frame(Y=Y, X=X, Z=Z, U=U, R=R)
  data = data_full
  # data[R==0, c("X.1","X.2")] = NA
  
  # output
  return(list(data = data, data_full = data_full))
}


#===========================================================================================#

# two functions for generating wish and iwish from peter hoff's book
rwish = function(n, nu0, S0) {
  sS0 = chol(S0)
  S = array(dim=c(dim(S0),n))
  for(i in 1:n) {
    Z = matrix(rnorm(nu0*dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i] = t(Z) %*% Z
  }
  S[,,1:n]
}

rinvwish = function(n,nu0,iS0) {
  sL0 = chol(iS0) 
  S = array(dim=c(dim(sL0),n))
  for(i in 1:n) {
    Z = matrix(rnorm(nu0 * dim(sL0)[1]), nu0, dim(iS0)[1]) %*% sL0  
    S[,,i] = solve(t(Z)%*%Z)
  }     
  S[,,1:n]
}
