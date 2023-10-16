// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_G(const arma::mat& v, const arma::sp_mat&  G, const double& prec=1.0e-8);
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& idx_q, const arma::vec& Xb_q, const arma::vec& off_q, const int& ns);
arma::sp_mat load_Q_hp(const arma::umat& from_to, const arma::vec& idx_q_r, const arma::vec& idx_q_m, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const arma::vec& off_q, const int& ns);
arma::mat load_L(const arma::vec& period_l, const arma::vec& cell_l, const arma::vec& idx_l, const arma::vec& fix_l, const arma::vec& Xb_l, const arma::vec& off_l, const int& ns, const int& np);

// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List mmpp_arma(const arma::vec& id, const  arma::vec& period, 
                     const arma::vec& dt, const arma::vec& cell, 
                     const int& ns, const int& np, 
                     const arma::mat& X_l, const arma::vec& off_l,
                     const arma::vec& fix_l, 
                     const arma::vec& period_l, const arma::vec& cell_l, 
                     const arma::vec& idx_l, 
                     const arma::vec& beta_l,
                     const arma::umat& from_to_q, const arma::mat& X_q_r, 
                     const arma::mat& X_q_m,
                     const arma::vec& off_q,
                     const arma::vec& idx_q_r, const arma::vec& idx_q_m,
                     const arma::vec& beta_q_r, const arma::vec& beta_q_m)
{
  int N = cell.size();
  
  arma::vec Xb_q_r = X_q_r * beta_q_r;
  arma::vec Xb_q_m = X_q_r * beta_q_m;
  arma::vec Xb_l = X_l * beta_l;

  double u = 0.0;
  arma::vec log_lik_v(N, fill::zeros);
  
  arma::sp_mat Q = load_Q_hp(from_to_q, idx_q_r, idx_q_m, Xb_q_r, Xb_q_m, off_q, ns);
  arma::sp_mat G(ns,ns);
  arma::mat L_mat = load_L(period_l, cell_l, idx_l, fix_l, Xb_l, off_l, ns, np);
  // Start forward loop
  arma::rowvec v(ns);
  arma::rowvec phi(ns);
  phi(cell(0)) = 1.0;
  arma::sp_mat P(ns,ns);
  arma::sp_mat L(ns,ns);
  
  // Start Forward alg loop (index = i)
  for(int i=1; i<N; i++){
    if(id(i)!=id(i-1)){
      phi = 0.0*phi; phi(cell(i)) = 1.0;
    } else{
      // Q = load_Q(from, to,Q_mat.col(period(i)), ns); #mod here for dynamic movement
      L.diag() = L_mat.col(period(i));
      G = (Q-L)*dt(i);
      v = phi_exp_G(phi, G);
      P = 0.0*P;
      if(R_finite(cell(i))){
        P(cell(i),cell(i)) = L(cell(i), cell(i));
        v = v * P; //elementwise
      }
      u = accu(v);
      log_lik_v(i) = log(u);
      phi = v/u;
    }
  } // end i
  
  double n2ll = -2*accu(log_lik_v);
  
  return Rcpp::List::create(
    // Rcpp::Named("log_lik_v") = log_lik_v,
    // Rcpp::Named("Q") = Q,
    // Rcpp::Named("G") = G,
    // Rcpp::Named("L_mat") = L_mat,
    // Rcpp::Named("L") = L,
    // Rcpp::Named("P") = P,
    // Rcpp::Named("v") = v,
    // Rcpp::Named("phi") = phi,
    // Rcpp::Named("N") = N,
    Rcpp::Named("n2ll") = n2ll
  );
  
}