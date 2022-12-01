// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <expQ2.h>

using namespace Rcpp;
using namespace expQ2;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat v_exp_QmLd(const arma::mat v, SEXP QmLd, double prec=1.0e-8) {
  arma::mat out;
  out = expQ2::v_exp_Q(v, QmLd, prec, false, true, false);
  return(out);
}

arma::mat phi_exp_QmLd(const arma::mat& v, const arma::sp_mat&  QmLd, double prec=1.0e-8) {
  arma::mat out = expQ2::sUnif_v_exp_Q(v,QmLd,prec,false,true);
  return out;
}


// [[Rcpp::export]]
arma::sp_mat load_Q(const arma::umat& from_to, const arma::uvec& idx_q, 
                    const arma::vec& qvals, const int& ns){
  int n = from_to.n_cols;
  arma::vec qext(n); 
  for(int i=0; i<n; i++){qext(i) = qvals(idx_q(i));}
  arma::sp_mat Q(from_to, qext, ns, ns);
  //arma::colvec ones(ns,fill::ones);
  arma::colvec row_sums = Q * ones(ns);
  Q.diag() = -1.0*row_sums;
  return Q;
}

// [[Rcpp::export]]
arma::mat load_L(const arma::uvec& period_l, const arma::uvec& cell_l, 
                 const arma::uvec& idx_l,
                 const arma::vec& fix_l, const arma::vec& l_vals, const int& ns, 
                 const int& np){
  arma::mat L_mat(ns,np);
  int n = period_l.size();
  for(int i=0; i<n; i++){
    if(!R_finite(fix_l(i))){L_mat(cell_l(i), period_l(i)) = l_vals(idx_l(i));
    } else{
      L_mat(cell_l(i), period_l(i)) = fix_l(i);
    }
  }
  return L_mat;
}



// Calculate likelihood ///////////////
// [[Rcpp::export]]
Rcpp::List mmpp_arma(const arma::vec& id, const  arma::uvec& period, 
                     const arma::vec& dt, const arma::uvec& cell, 
                     const int& ns, const int& np, 
                     const arma::mat& X_l, const arma::vec& fix_l, 
                     const arma::uvec& period_l, const arma::uvec& cell_l, 
                     const arma::uvec& idx_l, 
                     const arma::vec& beta_l,
                     const arma::umat& from_to_q, const arma::mat& X_q, 
                     const arma::uvec& idx_q, const arma::vec& beta_q)
{
  int N = cell.size();
  
  arma::vec q_vals = exp(X_q * beta_q);
  arma::vec l_vals = exp(X_l * beta_l);

  double u = 0.0;
  double log_lik = 0;
  arma::vec log_lik_v(N, fill::zeros);
  
  arma::sp_mat Q = load_Q(from_to_q, idx_q, q_vals, ns);
  arma::sp_mat QmLd(size(Q));
  arma::mat L_mat = load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np);
  // Start forward loop
  arma::mat v(1,ns);
  arma::mat phi(1,ns, fill::zeros);
  phi(0,cell(0)) = 1.0;
  arma::vec P(ns, fill::zeros);
  arma::vec L(ns);
  
  
  for(int i=1; i<N; i++){
    if(id(i)!=id(i-1)){
      phi = 0.0*phi; phi(0,cell(i)) = 1;
    } else{
      // Q = load_Q(from, to,Q_mat.col(period(i)), ns); #mod here for dynamic movement
      L = L_mat.col(period(i));
      QmLd = Q*dt(i); QmLd.diag() -= L * dt(i);
      v = phi_exp_QmLd(phi, QmLd);
      P = 0.0*P;
      if(R_finite(cell(i))){
        P(cell(i)) = L(cell(i));
        v = v % P; //elementwise
      } 
      u = accu(v);
      log_lik_v(i) = log(u);
      phi = v/u;
    }
  } // end i
  return Rcpp::List::create(
    Rcpp::Named("log_lik_v") = log_lik_v, 
    Rcpp::Named("Q") = Q,
    Rcpp::Named("L_mat") = L_mat,
    Rcpp::Named("L") = L,
    Rcpp::Named("P") = P,
    Rcpp::Named("v") = v,
    Rcpp::Named("phi") = phi,
    Rcpp::Named("N") = N
  );
  
}