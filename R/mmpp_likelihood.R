#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @author Devin S. Johnson
#' @import numDeriv
#' @export
mmpp_ll <- function(par, data_list,...){
from_to_q <- t(cbind(data_list$from, data_list$to))
beta_l <- par[1:ncol(data_list$X_l)]
beta_q <- par[(ncol(data_list$X_l)+1):(ncol(data_list$X_l)+ncol(data_list$X_q))]

mmpp_arma(
  data_list$id, 
  data_list$period, 
  data_list$dt, 
  data_list$cell, 
  data_list$ns, 
  data_list$np, 
  data_list$X_l, 
  data_list$fix_l, 
  data_list$period_l, 
  data_list$cell_l, 
  data_list$idx_l, 
  beta_l, 
  from_to_q, 
  data_list$X_q, 
  data_list$idx_q, 
  beta_q
  )$n2ll
}



