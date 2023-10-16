// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;



// [[Rcpp::export]]
arma::mat phi_exp_G(const arma::mat& v, const arma::sp_mat&  G, 
                    const double& prec=1.0e-8) {
  arma::mat out = expQ2::sv_exp_Q(v, G, prec, false, true);
  return out;
}


// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& idx_q, 
                    const arma::vec& Xb_q, const arma::vec& off_q, 
                    const int& ns){
  int n = from_to.n_cols;
  arma::vec qvals(n); 
  for(int i=0; i<n; i++){qvals(i) = exp(off_q(i) + Xb_q(idx_q(i)));}
  arma::sp_mat Q(from_to, qvals, ns, ns);
  //arma::colvec ones(ns,fill::ones);
  arma::colvec row_sums = Q * ones(ns);
  Q.diag() = -1.0*row_sums;
  return Q;
}



// [[Rcpp::export]]
arma::sp_mat load_Q_hp(const arma::umat& from_to, 
                       //const arma::vec& idx_q_r, const arma::vec& idx_q_m,
                       const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, 
                       const arma::vec& off_q, const int& ns){
  int n = from_to.n_cols;
  arma::sp_mat Qr(ns,ns);
  arma::vec Xb_m_vals(n);
  for(int i=0; i<ns; i++){Qr(i,i) = exp(off_q(i) + Xb_q_r(idx_q_r(i)));}
  for(int i=0; i<n; i++){Xb_m_vals(i) = Xb_q_m(idx_q_m(i));}
  arma::sp_mat Xb_m_mat(from_to, Xb_m_vals, ns, ns);
  arma::colvec x = max(Xb_m_mat, 1);
  
  Qm = normalise(Qm, 1, 1);
  arma::sp_mat Q = Qr * Qm;
  Q.diag() = -1.0*Qr.diag();
  return Q;
}


// [[Rcpp::export]]
arma::mat load_L(const arma::vec& period_l, const arma::vec& cell_l, 
                 //const arma::vec& idx_l,
                 const arma::vec& fix_l, const arma::vec& Xb_l, 
                 const arma::vec& off_l, const int& ns, 
                 const int& np){
  arma::mat L_mat(ns,np);
  int n = period_l.size();
  for(int i=0; i<n; i++){
    if(!R_finite(fix_l(i))){
      L_mat(cell_l(i), period_l(i)) = exp(off_l(i) + Xb_l(idx_l(i)));
    } else{
      L_mat(cell_l(i), period_l(i)) = fix_l(i);
    }
  }
  return L_mat;
}
