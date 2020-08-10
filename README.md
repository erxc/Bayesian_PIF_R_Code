# Bayesian PIF R Code

R code in this folder is for replicating simulation results in the Bayesian PIF estimation paper. 

## For the proposed estimation method:

1. To replicate simulation results with data generated from the bivariate normal model, use files simdata.R and chain_x_z_Rcpp.R. Parameter values can be modified in simdata.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp.R.

2. To replicate simulation results with data generated from the bivariate gamma model, use files simdata_mvgamma.R and chain_x_z_Rcpp.R. Parameter values can be modified in simdata_mvgamma.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp.R.

3. To replicate simulation results with data generated from the bivariate normal model and validation study participants selected as a biased sample, use files simdata_biased_sample.R and chain_x_z_Rcpp_biased_sample.R. Parameter values can be modified in simdata_biased_sample.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp_biased_sample.R.

## For the naive estimation method:

1. To replicate simulation results with data generated from the bivariate normal model, use files simdata.R and chain_x_z_Rcpp_naive.R. Parameter values can be modified in simdata.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp_naive.R.

2. To replicate simulation results with data generated from the bivariate gamma model, use files simdata_mvgamma.R and chain_x_z_Rcpp_naive.R. Parameter values can be modified in simdata_mvgamma.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp_naive.R.

3. To replicate simulation results with data generated from the bivariate normal model and validation study participants selected as a biased sample, use files simdata_biased_sample.R and chain_x_z_Rcpp_biased_sample_naive.R. Parameter values can be modified in simdata_biased_sample.R, while sizes of the main study and the validation study can be changes in chain_x_z_Rcpp_biased_sample_naive.R.
