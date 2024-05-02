#' @title Create model data necessary for MMPP simulation
#' @param ddl Design data list created from \code{\link[moveMMPP]{make_design_data}}.
#' @param model_parameters Model formula for the detection and movement portions of the MMPP model. 
#' @export
#' @author Devin S. Johnson
make_sim_dm <- function(ddl, model_parameters = mmpp_control()){
  
  if(max(ddl$lambda$period)>1) stop("Multiple time periods are not allowed in MMPP simulation at this time.")
  
  dml <- dm_lambda(model_parameters$lambda$form, ddl)
  dmq_r <- dm_q_r(model_parameters$q_r$form, ddl)
  dmq_m <- dm_q_m(model_parameters$q_m$form, ddl)
  
  par_map = list(
    beta_l = c(1:ncol(dml$X_l)),
    beta_q_r = c(1:ncol(dmq_r$X_q_r)) + ncol(dml$X_l)
  )
  if(ncol(dmq_m$X_q_m)!=0) par_map$beta_q_m = c(1:ncol(dmq_m$X_q_m)) + ncol(dmq_r$X_q_r) + ncol(dml$X_l)
  
  from_to <- ddl$q_m[,c("from_cellx","cellx")]
  fix_l <- ddl$lambda$fix
  
  return(list(X_l = dml$X_l, X_q_r=dmq_r$X_q_r, X_q_m=dmq_m$X_q_m, par_map=par_map, from_to=from_to, fix_l=fix_l, model_parameters=model_parameters))
}