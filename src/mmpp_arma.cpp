// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// function prototypes
arma::mat phi_exp_lnG(const arma::mat& phi, const arma::sp_mat&  lnG, const double& prec=1.0e-8);
arma::mat load_L(const arma::vec& period_l, const arma::vec& cell_l, const arma::vec& fix_l, const arma::vec& Xb_l, const int& ns, const int& np, const int& link_l=1, const double& a_l=1.0);
arma::sp_mat load_Q_mult(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const int& link_m=1, const double& a_r=1.0, const double& a_m=1.0, const bool& norm=true);
arma::sp_mat load_Q_add(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const int& link_r=1, const int& link_m=1, const double& a_r=1.0, const double& a_m=1.0);
arma::sp_mat load_Q_sde(const arma::umat& from_to, const arma::vec& Xb_q_r, const arma::vec& Xb_q_m, const int& ns, const double& k,const double& a_r=1.0);


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
                     const double& eq_prec = 1.0e-8,
                     const int& link_l = 1,  
                     const int& link_r = 1,
                     const int& link_m = 1,
                     const int& struc = 1,
                     const double& a_l = 1.0,
                     const double& a_r = 1.0, 
                     const double& a_m = 1.0,
                     const bool& norm=true)
{
  int N = cell.size();

  double u = 0.0;
  arma::vec log_lik_v(N, fill::zeros);
  
  arma::sp_mat Q;
  if(struc==1){
    Q = load_Q_mult(from_to, Xb_q_r, Xb_q_m, ns, link_r, link_m, a_r, a_m, norm);
  } else { //}if(struc==2){
    Q = load_Q_add(from_to, Xb_q_r, Xb_q_m, ns, link_r, link_m, a_r, a_m);
  } 
  // else if(struc==3){
  //   Q = load_Q_sde(from_to, Xb_q_r, Xb_q_m, ns, k, a_r);
  // }
  
  arma::sp_mat lnG(ns,ns);
  arma::mat L_mat = load_L(period_l, cell_l, fix_l, Xb_l, ns, np, link_l, a_l);
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