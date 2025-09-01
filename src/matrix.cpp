// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::mat solve_cpp(arma::mat M,bool symmetric){
  arma::mat M_inv(arma::size(M));
  if(symmetric){

    M_inv= arma::inv_sympd(M);
  }
  else{
    M_inv= arma::inv(M);
  }
  return M_inv;
}


// [[Rcpp::export]]
arma::sp_mat spsolve_cpp(arma::sp_mat M){
  arma::sp_mat M_inv(arma::size(M));
  arma::mat Id(arma::size(M));
  M_inv=arma::spsolve(M,Id);

  return M_inv;
}

// [[Rcpp::export]]
arma::mat kronecker_cpp(arma::mat M,arma::mat Q){
  return arma::kron(M,Q);
}

// [[Rcpp::export]]
arma::sp_mat spkronecker_cpp(arma::sp_mat M,arma::sp_mat Q){
  return arma::kron(M,Q);
}