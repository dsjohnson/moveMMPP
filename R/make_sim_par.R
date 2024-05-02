#' @title Create Q and L matrices for MMPP simulation
#' @param sim_dm A list of design matrices resulting from a call to \code{\link[moveMMPP]{make_sim_dm}}
#' @param par A list corresponding to the coefficients for lambda, Q_m, and Q_r. It must match the 
#' lengths in \code{sim_dm$par_map}.
#' @export
make_sim_par <- function(sim_dm, par){
  par_map = sim_dm$par_map
  model_parameters <- sim_dm$model_parameters
  
  if(length(par$beta_l) != length(par_map$beta_l)) stop("Length of beta_l vector is not correct!")
  if(length(par$beta_q_r) != length(par_map$beta_q_r)) stop("Length of beta_q_r vector is not correct!")
  if(length(par$beta_q_m) != length(par_map$beta_q_m)) stop("Length of beta_q_m vector is not correct!")
  
  from_to <- t(as.matrix(sim_dm$from_to))-1
  ns <- nrow(sim_dm$X_l)
  Xb_l <- sim_dm$X_l %*% par$beta_l
  Xb_q_r <- sim_dm$X_q_r %*% par$beta_q_r
  Xb_q_m <-sim_dm$X_q_m %*% par$beta_q_m
  link_r <- model_parameters$q_r$link
  a_r <-  model_parameters$q_r$a
  link_m <- model_parameters$q_m$link
  a_m <-  model_parameters$q_m$a
  norm <- model_parameters$norm
  
  
  if(model_parameters$struc=="mult"){
    Q <- load_Q_mult(from_to, Xb_q_r, Xb_q_m, ns, 
                     link_r = which(link_r==c("soft_plus", "log")), 
                     link_m = which(link_m==c("soft_plus", "log")), 
                     a_r, a_m, norm)
  } else {
    Q <- load_Q_add(from_to, Xb_q_r, Xb_q_m, ns,  
                    link_r = which(link_r==c("soft_plus", "log")), 
                    link_m = which(link_m==c("soft_plus", "log")), 
                    a_r, a_m) 
  }
  L <- load_L(period_l=rep(0,ns), cell_l=c(1:ns)-1, fix_l=sim_dm$fix_l, Xb_l=Xb_l, ns=ns, np=1)
  
  return(list(Q=Q, L=L))
  
}