# Bayesian PIF `R` Code

The `R` code in this folder is used to replicate simulation results in the Bayesian PIF estimation paper. The exposure-outcome relationship is assumed to be a logistic regression, with the default prevalence being -2 or -4. The reclassification process is assumed to follow a multivariate normal model. For different simulation scenarios mentioned in the paper, the default prevalence needs to be changed in both the data-generating file and the MCMC chain file. The reclassification process parameters can be altered in the data generation file to achieve different levels of measurement errors. Sizes of the main study and validation study can be altered in MCMC chain files.

## For the proposed estimation method:

1. To replicate simulation results with data generated from the multivariate normal model, use files `simdata.R` and `chain_x_z_Rcpp.R`. Parameter values can be modified in `simdata.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp.R`.

2. To replicate simulation results with data generated from the multivariate gamma model, use files `simdata_mvgamma.R` and `chain_x_z_Rcpp.R`. Parameter values can be modified in `simdata_mvgamma.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp.R`.

3. To replicate simulation results with data generated from the multivariate normal model and validation study participants selected as a biased sample, use files `simdata_biased_sample.R` and `chain_x_z_Rcpp_biased_sample.R`. Parameter values can be modified in `simdata_biased_sample.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp_biased_sample.R`.

## For the naive estimation method:

1. To replicate simulation results with data generated from the multivariate normal model, use files `simdata.R` and `chain_x_z_Rcpp_naive.R`. Parameter values can be modified in `simdata.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp_naive.R`.

2. To replicate simulation results with data generated from the multivariate gamma model, use files `simdata_mvgamma.R` and `chain_x_z_Rcpp_naive.R`. Parameter values can be modified in `simdata_mvgamma.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp_naive.R`.

3. To replicate simulation results with data generated from the multivariate normal model and validation study participants selected as a biased sample, use files `simdata_biased_sample.R` and `chain_x_z_Rcpp_biased_sample_naive.R`. Parameter values can be modified in `simdata_biased_sample.R`, while the sizes of the main study and the validation study can be changed in `chain_x_z_Rcpp_biased_sample_naive.R`.
