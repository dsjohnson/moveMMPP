// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <iostream>
using namespace arma;

// [[Rcpp::export]]
double sum_by_iterator(const arma::sp_mat& x) {
  double result = 0;
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    result += *i * (i.row() + 1);
  }
  return result;
}

// [[Rcpp::export]]
arma::sp_mat sweep_max(const arma::sp_mat& X) {
  arma::sp_mat out = X;
  arma::vec mx;
  mx = arma::max(X, 1);
  for (arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i) {
    out(i.row(), i.col()) = *i - mx(i.row());
  }
  return out;
}

// [[Rcpp::export]]
arma::sp_mat sp_exp(const arma::sp_mat& X, const bool& row_sweep = true) {
  arma::sp_mat out = X;
  if(row_sweep){
    arma::vec mx;
    mx = arma::max(X, 1);
    for (arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i) {
      out(i.row(), i.col()) = exp(*i - mx(i.row()));
    }
  } else{
    for (arma::sp_mat::const_iterator i = X.begin(); i != X.end(); ++i) {
      out(i.row(), i.col()) = exp(*i);
    }
  } 
  return out;
}


// [[Rcpp::export]]
arma::sp_mat sp_norm(const arma::sp_mat& X) {
  arma::sp_mat out = X;
  out =  normalise(X, 1, 1);
  out.diag().ones();
  return out;
}


// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to,
                    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                    const int& ns, const bool& row_sweep=true) {
  arma::sp_mat Qr(ns,ns);
  Qr.diag() = trunc_exp(Xb_q_r);
  arma::sp_mat Xb_m_mat(from_to, Xb_q_m, ns, ns);
  // Xb_m_mat = sp_sweep_max( Xb_m_mat);
  arma::sp_mat Qm = sp_exp(Xb_m_mat, row_sweep);
  Qm = normalise(Qm, 1, 1);
  Qm.diag().ones();
  Qm.diag() = -1*Qm.diag();
  arma::sp_mat Q = Qr * Qm;
  return Q;
}

// [[Rcpp::export]]
arma::sp_mat load_Q2(const arma::umat& from_to,
                    const arma::vec& Q_r, const arma::vec& Xb_q_m,
                    const int& ns, const bool& row_sweep=true) {
  arma::sp_mat Qr(ns,ns);
  Qr.diag() = Q_r;
  arma::sp_mat Xb_m_mat(from_to, Xb_q_m, ns, ns);
  // Xb_m_mat = sp_sweep_max( Xb_m_mat);
  arma::sp_mat Qm = sp_exp(Xb_m_mat, row_sweep);
  Qm = normalise(Qm, 1, 1);
  Qm.diag().ones();
  Qm.diag() = -1*Qm.diag();
  arma::sp_mat Q = Qr * Qm;
  return Q;
}

