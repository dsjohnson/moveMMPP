// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const bool& row_sweep=true);
arma::mat load_L(const arma::vec& period_l, const arma::vec& cell_l, const arma::vec& fix_l, const arma::vec& Xb_l, const int& ns, const int& np);

// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List mmpp_arma(const arma::vec& id, const  arma::vec& period, 
                     const arma::vec& dt, const arma::vec& cell, 
                     const int& ns, const int& np, 
                     const arma::vec& Xb_l, 
                     const arma::vec& fix_l, 
                     const arma::vec& period_l, const arma::vec& cell_l, 
                     const arma::umat& from_to, 
                     const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                     const bool& row_sweep=true)
{
  int N = cell.size();
  
  //arma::vec Xb_q_r = X_q_r * beta_q_r;
  //arma::vec Xb_q_m = X_q_m * beta_q_m;
  //arma::vec Xb_l = X_l * beta_l;
  
  // arma::vec q_vals = exp();
  // arma::vec l_vals = exp();

  double u = 0.0;
  arma::vec log_lik_v(N, fill::zeros);
  
  arma::sp_mat Q = load_Q(from_to, Xb_q_r, Xb_q_m, ns, row_sweep);
  arma::sp_mat lnG(ns,ns);
  arma::mat L_mat = load_L(period_l, cell_l, fix_l, Xb_l, ns, np);
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
      lnG = (Q-L)*dt(i);
      v = phi_exp_lnG(phi, lnG);
      P.zeros(); 
      if(R_finite(cell(i))){
        P(cell(i),cell(i)) = L(cell(i), cell(i));
        v = v * P; 
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