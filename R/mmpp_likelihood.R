#' @title Evaluate movement MMPP log-likelihood
#' @param par Parameter vector
#' @param data_list List of required data objects to evaluate likelihood
#' @param debug For developers only, leave in the default setting.
#' @param ... Extra wiggle room for ignored arguments.
#' @importFrom stats aggregate
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

  
  # if(is.character(data_list$delta) && data_list$delta=="stationary"){
  #   delta <- get_lim_ud(list(par = par, data_list = data_list))
  #   delta <- delta$ud
  #   delta <- delta/sum(delta)
  # } else{
  #   delta <- data_list$delta
  # }
  
  
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
    Xb_q_m,
    # delta = matrix(delta, nrow=1),
    eq_prec = data_list$eq_prec,
    link_r = which(data_list$link_r==c("soft_plus", "log", "logit")),
    link_m = which(data_list$link_m==c("soft_plus", "log", "logit")),
    struc = which(data_list$struc==c("mult", "add", "sde")),
    a_r = data_list$a_r,
    a_m = data_list$a_m,
    norm = data_list$norm
  )$n2ll
}



