#' @title Fit MMPP movement model to resight data
#' @param data A processed data frame produced by the function \code{\link{process_data}}
#' @param ddl A design data list produced by the function \code{\link{make_design_data}}.
#' @param model_parameters Model formula for the detection and movement portions
#' of the MMPP model. Must be of the form \code{list(lambda=~form, q=~form)}.
#' @param fit_model Logical. Should the model be fit, or just return the TMB object? 
#' Defaults to \code{fit_model = TRUE}.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q=c())}.
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{opm}} from the \code{\link[optimx]{optimx-package}}.
#' @author Devin S. Johnson
#' @import optimx TMB dplyr
#' @export
fit_mmpp <- function(data, ddl, 
                     model_parameters=list(lambda=~1, q=~1),
                     fit_model=TRUE,
                     start=NULL, ...
                     ){
  cell <- cellx <- NULL
  
  cell_idx_df <- select(ddl$lambda, cell, cellx) %>% distinct()
  data <- data %>% left_join(cell_idx_df, by="cell")
  
  dml_list <- dm_lambda(model_parameters$lambda, ddl)
  dmq_list <- dm_q(model_parameters$q, ddl)
  
  tmb_data <- list(
    N = as.integer(nrow(data)),
    ns = as.integer(length(unique(ddl$lambda$cell))),
    np = as.integer(max(ddl$quad_pts$period)),
    #detections
    id = data$idx-1,
    period = as.integer(data$period-1),
    dt = data$delta,
    cell = as.integer(data$cellx-1),
    # lambda
    X_l = dml_list$X_l,
    fix_l = dml_list$idx_l$fix,
    period_l = as.integer(dml_list$idx_l$period-1),
    cell_l = as.integer(dml_list$idx_l$cell-1),
    idx_l = as.integer(dml_list$idx_l$idx_l-1),
    # Q
    from_q = as.integer(dmq_list$idx_q$from_cellx-1),
    to_q = as.integer(dmq_list$idx_q$to_cellx-1),
    X_q = dmq_list$X_q,
    idx_q = as.integer(dmq_list$idx_q$idx_q-1)
  )
  
  if(!is.null(start)){
  tmb_par <- list(
    beta_l=rep(0,ncol(dml_list$X_l)) , 
    beta_q=rep(0,ncol(dmq_list$X_q))
    )
  } else{
    tmb_par=start
  }
  
  foo <- MakeADFun(
    data=append(list(model="mmpp"), tmb_data),
    parameters=tmb_par,
    #random=c(),
    DLL="moveMMPP_TMBExports"
  )
  
opt <- opm(f$par,f$fn,f$gr,f$he,hessian=FALSE,...)
  
}