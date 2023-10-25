#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param ... Extra wiggle room for ignored arguments.
#' @author Devin S. Johnson
#' @export
mmpp_ll <- function(par, data_list,...){
from_to <- t(cbind(data_list$from, data_list$to))
par_map <- data_list$par_map
beta_l <- par[par_map$beta_l]
beta_q_r <- par[par_map$beta_q_r]
beta_q_m <- par[par_map$beta_q_m]

Xb_l <- data_list$X_l %*% beta_l
Xb_q_r <- data_list$X_q_r %*% beta_q_r
Xb_q_m <- data_list$X_q_m %*% beta_q_m

if(all(Xb_q_m == 0) & ncol(X_q_m)==0) X_q_m = X_q_m + 1

mmpp_arma(
  data_list$id, 
  data_list$period, 
  data_list$dt, 
  data_list$cell, 
  data_list$ns, 
  data_list$np, 
  Xb_l, 
  # data_list$off_l,
  data_list$fix_l, 
  data_list$period_l, 
  data_list$cell_l, 
  # data_list$idx_l, 
  from_to, 
  Xb_q_r, 
  Xb_q_m #,
  # data_list$off_q,
  # data_list$idx_q, 
  )$n2ll
}



