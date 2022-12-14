// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// phi_exp_G
arma::mat phi_exp_G(const arma::mat& v, const arma::sp_mat& G, const double& prec);
RcppExport SEXP _moveMMPP_phi_exp_G(SEXP vSEXP, SEXP GSEXP, SEXP precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type G(GSEXP);
    Rcpp::traits::input_parameter< const double& >::type prec(precSEXP);
    rcpp_result_gen = Rcpp::wrap(phi_exp_G(v, G, prec));
    return rcpp_result_gen;
END_RCPP
}
// load_Q
arma::sp_mat load_Q(const arma::umat& from_to, const arma::vec& idx_q, const arma::vec& Xb_q, const arma::vec& off_q, const int& ns);
RcppExport SEXP _moveMMPP_load_Q(SEXP from_toSEXP, SEXP idx_qSEXP, SEXP Xb_qSEXP, SEXP off_qSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type from_to(from_toSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx_q(idx_qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xb_q(Xb_qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type off_q(off_qSEXP);
    Rcpp::traits::input_parameter< const int& >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(load_Q(from_to, idx_q, Xb_q, off_q, ns));
    return rcpp_result_gen;
END_RCPP
}
// load_L
arma::mat load_L(const arma::vec& period_l, const arma::vec& cell_l, const arma::vec& idx_l, const arma::vec& fix_l, const arma::vec& Xb_l, const arma::vec& off_l, const int& ns, const int& np);
RcppExport SEXP _moveMMPP_load_L(SEXP period_lSEXP, SEXP cell_lSEXP, SEXP idx_lSEXP, SEXP fix_lSEXP, SEXP Xb_lSEXP, SEXP off_lSEXP, SEXP nsSEXP, SEXP npSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type period_l(period_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cell_l(cell_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx_l(idx_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fix_l(fix_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Xb_l(Xb_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type off_l(off_lSEXP);
    Rcpp::traits::input_parameter< const int& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int& >::type np(npSEXP);
    rcpp_result_gen = Rcpp::wrap(load_L(period_l, cell_l, idx_l, fix_l, Xb_l, off_l, ns, np));
    return rcpp_result_gen;
END_RCPP
}
// mmpp_arma
Rcpp::List mmpp_arma(const arma::vec& id, const arma::vec& period, const arma::vec& dt, const arma::vec& cell, const int& ns, const int& np, const arma::mat& X_l, const arma::vec& off_l, const arma::vec& fix_l, const arma::vec& period_l, const arma::vec& cell_l, const arma::vec& idx_l, const arma::vec& beta_l, const arma::umat& from_to_q, const arma::mat& X_q, const arma::vec& off_q, const arma::vec& idx_q, const arma::vec& beta_q);
RcppExport SEXP _moveMMPP_mmpp_arma(SEXP idSEXP, SEXP periodSEXP, SEXP dtSEXP, SEXP cellSEXP, SEXP nsSEXP, SEXP npSEXP, SEXP X_lSEXP, SEXP off_lSEXP, SEXP fix_lSEXP, SEXP period_lSEXP, SEXP cell_lSEXP, SEXP idx_lSEXP, SEXP beta_lSEXP, SEXP from_to_qSEXP, SEXP X_qSEXP, SEXP off_qSEXP, SEXP idx_qSEXP, SEXP beta_qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type period(periodSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cell(cellSEXP);
    Rcpp::traits::input_parameter< const int& >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< const int& >::type np(npSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_l(X_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type off_l(off_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fix_l(fix_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type period_l(period_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cell_l(cell_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx_l(idx_lSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_l(beta_lSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type from_to_q(from_to_qSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X_q(X_qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type off_q(off_qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type idx_q(idx_qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta_q(beta_qSEXP);
    rcpp_result_gen = Rcpp::wrap(mmpp_arma(id, period, dt, cell, ns, np, X_l, off_l, fix_l, period_l, cell_l, idx_l, beta_l, from_to_q, X_q, off_q, idx_q, beta_q));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_moveMMPP_phi_exp_G", (DL_FUNC) &_moveMMPP_phi_exp_G, 3},
    {"_moveMMPP_load_Q", (DL_FUNC) &_moveMMPP_load_Q, 5},
    {"_moveMMPP_load_L", (DL_FUNC) &_moveMMPP_load_L, 8},
    {"_moveMMPP_mmpp_arma", (DL_FUNC) &_moveMMPP_mmpp_arma, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_moveMMPP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
