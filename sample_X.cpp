
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "f_mvnorm.h"
#include "f_madd.h"
#include "f_mmult.h"
#include "f_inv.h"

using namespace Rcpp;

// NumericMatrix

// [[Rcpp::export]]
NumericMatrix sample_X(int n_m_cpp, int d_x_cpp, NumericVector beta_x_cpp, NumericMatrix c_beta_x_cpp,
                       NumericVector w_cpp, NumericVector eta_nx_cpp, NumericMatrix zeta_cpp, 
                       NumericMatrix invSigma_cpp, NumericVector Y_cpp) {
  NumericMatrix X_mis_cpp(n_m_cpp, d_x_cpp);

  for (int i1=0; i1<n_m_cpp; ++i1) {
    
    NumericMatrix V_i(d_x_cpp, d_x_cpp);
    NumericMatrix V_i_tmp(d_x_cpp, d_x_cpp);
    NumericVector mu_i(d_x_cpp);
    NumericMatrix mu_i_tmp(d_x_cpp, 1);
    NumericMatrix t_zeta_k(d_x_cpp, 1);
    NumericMatrix inv_Sig_zeta(d_x_cpp, 1);
    NumericMatrix y_w_eta_beta(d_x_cpp, 1);
    NumericMatrix x_new(d_x_cpp, 1);
    
    for (int k1=0; k1<d_x_cpp; ++k1) {
      for (int k2=0; k2<d_x_cpp; ++k2) {
        V_i_tmp(k1,k2) = w_cpp[i1]*c_beta_x_cpp(k1,k2);
      }
    }
    
    V_i = f_inv(f_madd(V_i_tmp, invSigma_cpp));
    
    for (int j1=0; j1<d_x_cpp; ++j1) {
      t_zeta_k(j1,0) = zeta_cpp(i1,j1);
    }
    
    inv_Sig_zeta = f_mmult(invSigma_cpp, t_zeta_k);
    
    for (int j2=0; j2<d_x_cpp; ++j2) {
      y_w_eta_beta(j2,0) = (Y_cpp[i1]-0.5-w_cpp[i1]*eta_nx_cpp[i1]) * beta_x_cpp[j2];
    }
    
    mu_i_tmp = f_mmult(V_i, f_madd(y_w_eta_beta, inv_Sig_zeta));
    
    for (int j3=0; j3<d_x_cpp; ++j3) {
      mu_i[j3] = mu_i_tmp(j3,0);
    }
    
    x_new = f_mvnorm(1, mu_i, V_i);

    for (int j4=0; j4<d_x_cpp; ++j4) {
      X_mis_cpp(i1,j4) = x_new(j4,0);
    }
  }
  
  return(X_mis_cpp);
  // return List::create(Rcpp::Named("Gam_x") = Gam_x,
  //                     Rcpp::Named("Gam_x_t") = Gam_x_t,
  //                     Rcpp::Named("V_i_tmp") = V_i_tmp,
  //                     Rcpp::Named("beta_o") = beta_o,
  //                     Rcpp::Named("beta_o_w") = TEMP);
}