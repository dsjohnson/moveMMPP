#include <RcppArmadillo.h>
#include <expQ2.h>
#include "../inst/include/moveMMPP.h"

using namespace Rcpp;
using namespace arma;

// Declaration of internal function from mmpp_arma.cpp
Rcpp::List mmpp_arma(const arma::vec& id, const arma::vec& period, 
                     const arma::vec& dt, const arma::vec& cell, 
                     const int& ns, const int& np, 
                     const arma::vec& Xb_l, const arma::vec& fix_l, 
                     const arma::vec& period_l, const arma::vec& cell_l, 
                     const arma::umat& from_to, 
                     const arma::vec& Xb_q_r, const arma::vec& Xb_q_m,
                     const double& eq_prec,
                     const int& link_l, const int& link_r, const int& link_m, const int& struc,
                     const double& a_l, const double& a_r, const double& a_m, const bool& norm);

extern "C" {
  
  double c_mmpp_loglik(
      double* id_ptr, int n_id,
      double* period_ptr, int n_period,
      double* dt_ptr, int n_dt,
      double* cell_ptr, int n_cell,
      int ns, int np,
      double* Xb_l_ptr, int n_Xb_l,
      double* fix_l_ptr, int n_fix_l,
      double* period_l_ptr, int n_period_l,
      double* cell_l_ptr, int n_cell_l,
      double* from_to_ptr, int n_from_to_rows, int n_from_to_cols,
      double* Xb_q_r_ptr, int n_Xb_q_r,
      double* Xb_q_m_ptr, int n_Xb_q_m,
      double eq_prec,
      int link_l, int link_r, int link_m, int struc,
      double a_l, double a_r, double a_m, int norm_int,
      int log_p
  ) {
    // Reconstitute vectors via zero-copy constructors
    arma::vec id(id_ptr, n_id, false, true);
    arma::vec period(period_ptr, n_period, false, true);
    arma::vec dt(dt_ptr, n_dt, false, true);
    arma::vec cell(cell_ptr, n_cell, false, true);
    
    arma::vec Xb_l(Xb_l_ptr, n_Xb_l, false, true);
    arma::vec fix_l(fix_l_ptr, n_fix_l, false, true);
    arma::vec period_l(period_l_ptr, n_period_l, false, true);
    arma::vec cell_l(cell_l_ptr, n_cell_l, false, true);
    
    arma::vec Xb_q_r(Xb_q_r_ptr, n_Xb_q_r, false, true);
    arma::vec Xb_q_m(Xb_q_m_ptr, n_Xb_q_m, false, true);
    
    arma::mat from_to_dbl(from_to_ptr, n_from_to_rows, n_from_to_cols, false, true);
    arma::umat from_to = arma::conv_to<arma::umat>::from(from_to_dbl);
    
    bool norm = (norm_int != 0);
    
    // Compute likelihood
    Rcpp::List res = mmpp_arma(
      id, period, dt, cell, ns, np,
      Xb_l, fix_l, period_l, cell_l,
      from_to, Xb_q_r, Xb_q_m,
      eq_prec, link_l, link_r, link_m, struc,
      a_l, a_r, a_m, norm
    );
    
    double n2ll = Rcpp::as<double>(res["n2ll"]);
    double loglik = -0.5 * n2ll;
    
    return log_p ? loglik : std::exp(loglik);
  }
  
}