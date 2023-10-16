#' @title Evaluate movement MMPP log-likelihood with Hewitt movement rate form
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param ... Extra wiggle room for ignored arguments.
#' @author Devin S. Johnson
#' @export
mmpp_ll_hp <- function(par, data_list,...){
from_to_q <- t(cbind(data_list$from, data_list$to))
beta_l <- par[1:ncol(data_list$X_l)]
beta_q_r <- par[(ncol(data_list$X_l)+1):(ncol(data_list$X_l)+ncol(data_list$X_q))]
beta_q_m <- par[(ncol(data_list$X_l)+1):(ncol(data_list$X_l)+ncol(data_list$X_q))]

mmpp_arma(
  data_list$id, 
  data_list$period, 
  data_list$dt, 
  data_list$cell, 
  data_list$ns, 
  data_list$np, 
  data_list$X_l, 
  data_list$off_l,
  data_list$fix_l, 
  data_list$period_l, 
  data_list$cell_l, 
  data_list$idx_l, 
  beta_l, 
  from_to_q, 
  data_list$X_q, 
  data_list$off_q,
  data_list$idx_q, 
  beta_q
  )$n2ll
}



