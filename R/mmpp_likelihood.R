#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @author Devin S. Johnson
#' @export
mmpp_ll <- function(par, data_list, debug=0, ...){
  if(debug>0) browser()
  from_to <- t(cbind(data_list$from, data_list$to))
  par_map <- data_list$par_map
  beta_l <- par[par_map$beta_l]
  beta_q_r <- par[par_map$beta_q_r]
  beta_q_m <- par[par_map$beta_q_m]
  
  Xb_l <- data_list$X_l %*% beta_l
  Xb_q_r <- data_list$X_q_r %*% beta_q_r
  Xb_q_m <- data_list$X_q_m %*% beta_q_m
  # e <- rnorm(length(Xb_q_m))
  # Xb_q_m = Xb_q_m + e
  
  #if(ncol(data_list$X_q_m)==0) Xb_q_m = Xb_q_m + 1
  
  mx <- (aggregate(Xb_q_m, list(data_list$from), max)[,2])[data_list$from+1]
  Xb_q_m <- Xb_q_m-mx
  
  # L_mat <- load_L(data_list$period_l, data_list$cell_l, data_list$fix_l, Xb_l, data_list$ns,data_list$np)
  # qqq <- load_Q(from_to, Xb_q_r, Xb_q_m, data_list$ns, norm=TRUE)
  
  mmpp_arma(
    data_list$id, 
    data_list$period, 
    data_list$dt, 
    data_list$cell, 
    data_list$ns, 
    data_list$np, 
    Xb_l, 
    data_list$fix_l, 
    data_list$period_l, 
    data_list$cell_l, 
    from_to, 
    Xb_q_r, 
    Xb_q_m
  )$n2ll
}



