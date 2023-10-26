#define arma_64bit_word 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[rcpp::plugins(cpp11)]] 
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;



// [[Rcpp::export]]
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(phi, lnG, prec, false, true);
  return out;
}


// // [[Rcpp::export]]
// arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& idx_q, 
//                     const arma::vec& Xb_q, const arma::vec& off_q, 
//                     const int& ns){
//   int n = from_to.n_cols;
//   arma::vec qvals(n); 
//   for(int i=0; i<n; i++){qvals(i) = exp(off_q(i) + Xb_q(idx_q(i)));}
//   arma::sp_mat Q(from_to, qvals, ns, ns);
//   //arma::colvec ones(ns,fill::ones);
//   arma::colvec row_sums = Q * ones(ns);
//   Q.diag() = -1.0*row_sums;
//   return Q;
// }


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



// // [[Rcpp::export]]
// arma::sp_mat load_Q(const arma::umat& from_to,
//                     const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
//                     const int& ns, const bool& row_sweep=true) {
//   arma::sp_mat Qr(ns,ns);
//   Qr.diag() = trunc_exp(Xb_q_r);
//   arma::sp_mat Xb_m_mat(from_to, Xb_q_m, ns, ns);
//   arma::sp_mat Qm = sp_exp(Xb_m_mat, row_sweep);
//   Qm = normalise(Qm, 1, 1);
//   Qm.diag().ones();
//   Qm.diag() = -1*Qm.diag();
//   arma::sp_mat Q = Qr * Qm;
//   return Q;
// }

// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to,
                    const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                    const int& ns, const bool& norm=true) {
  arma::sp_mat Qr(ns,ns);
  Qr.diag() = trunc_exp(Xb_q_r);
  arma::sp_mat Qm(from_to, trunc_exp(Xb_q_m), ns, ns);
  if(norm){
    Qm = normalise(Qm, 1, 1);
  }
  Qm.diag().ones();
  Qm.diag() = -1*Qm.diag();
  arma::sp_mat Q = Qr * Qm;
  return Q;
}

// [[Rcpp::export]]
arma::mat load_L(
    const arma::vec& period_l,
    const arma::vec& cell_l,
    const arma::vec& fix_l,
    const arma::vec& Xb_l,
    const int& ns, const int& np
){
  arma::mat L_mat(ns,np);
  int n = period_l.size();
  for(int i=0; i<n; i++){
    if(!R_finite(fix_l(i))){
      L_mat(cell_l(i), period_l(i)) = trunc_exp(Xb_l(i));
    } else{
      L_mat(cell_l(i), period_l(i)) = fix_l(i);
    }
  }
  return L_mat;
}