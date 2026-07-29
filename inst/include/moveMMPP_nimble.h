#ifndef MOVEMMPP_NIMBLE_H
#define MOVEMMPP_NIMBLE_H

#ifdef __cplusplus
extern "C" {
#endif
  
  // Plain C signature for NIMBLE (no Rcpp/SEXP dependency)
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
  );
  
#ifdef __cplusplus
}
#endif

#endif // MOVEMMPP_NIMBLE_H