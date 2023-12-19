#' @title Fit MMPP movement model to resight data
#' @param data A processed data frame produced by the function \code{\link{process_data}}
#' @param ddl A design data list produced by the function \code{\link{make_design_data}}.
#' @param model_parameters Model formula for the detection and movement portions
#' of the MMPP model. 
#' @param pen_fun An optional penalty function. Should be on the scale of a log-prior distribution.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q_r=c(), beta_q_r=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param fit Logical. Should the likelihood be optimized?
#' @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#' package developers. 
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @details
#' Two model forms are available \code{list(lambda=list(form=~1, offset=NULL), q=list(form=~1, offset=~log(1.num_neigh)))}. For the 
#' \code{q} model you must use \code{offset=0} to not have one. If it is left off, \code{offset=~log(1.num_neigh)))} 
#' will be used. To use the movement rate form of Hewitt et al. (2023), one must use, e.g., 
#' `q = list(res_form = ~from_var, mov_form=~to_var, offset=...)`, where `from_var` is a variable
#' @references Hewitt, J., Gelfand, A. E., & Schick, R. S. (2023). Time-discretization approximation enriches continuous-time discrete-space models for animal movement. The Annals of Applied Statistics, 17:740-760.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @importFrom stats ppois
#' @export
fit_mmpp_dir <- function(data, ddl, 
                     model_parameters = list(lambda = ~1, q_r = ~1, q_m = ~1
                     ), pen_fun = NULL,
                     hessian=TRUE, start=NULL, method="nlminb", fit=TRUE, 
                     debug=0, ...){
  
  cell <- cellx <- fix <- NULL
  
  if(debug==1) browser()
  
  cell_idx_df <- select(ddl$q_r, cell, cellx) %>% distinct()
  data <- data %>% left_join(cell_idx_df, by="cell")
  
  dml <- dm_lambda(model_parameters$lambda, ddl)
  dmq_r <- dm_q_r(model_parameters$q_r, ddl)
  dmq_m <- dm_q_m(model_parameters$q_m, ddl)
  
  data$period <- ifelse(is.na(data$cell), data$period-1, data$period)
  
  par_map = list(
    beta_l = c(1:ncol(dml$X_l)),
    beta_q_r = c(1:ncol(dmq_r$X_q_r)) + ncol(dml$X_l)
  )
  if(ncol(dmq_m$X_q_m)!=0) par_map$beta_q_m = c(1:ncol(dmq_m$X_q_m)) + ncol(dmq_r$X_q_r) + ncol(dml$X_l)

  
  data_list <- list(
    N = as.integer(nrow(data)),
    ns = as.integer(length(unique(ddl$q_r$cellx))),
    np = as.integer(max(ddl$quad_pts$period)),
    # detection
    id = data$idx-1,
    period = as.integer(data$period-1),
    dt = data$delta,
    cell = as.integer(data$cellx-1),
    ### lambda
    X_l = dml$X_l,
    # off_l = dml_list$off_l,
    fix_l = ddl$lambda$fix,
    period_l = as.integer(ddl$lambda$period-1),
    cell_l = as.integer(ddl$lambda$cellx-1),
    # idx_l = as.integer(dml_list$idx_l$idx_l-1),
    ### Q
    from = as.integer(ddl$q_m$from_cellx-1),
    to = as.integer(ddl$q_m$cellx-1),
    X_q_r = dmq_r$X_q_r,
    X_q_m = dmq_m$X_q_m,
    # off_q = dmq_list$off_q,
    # idx_q = as.integer(dmq_list$idx_q$idx_q-1)
    par_map = par_map
  )
  
  if(is.null(start)){
    par_list <- list(
      beta_l = rep(0,ncol(dml$X_l)), 
      beta_q_r = rep(0, ncol(dmq_r$X_q_r)),
      beta_q_m = rep(0, ncol(dmq_m$X_q_m))
    )
  } else{
    par_list=start
  }
  
  start <- c(par_list$beta_l, par_list$beta_q_r, par_list$beta_q_m)
  
  if(is.null(pen_fun)){
    obj_fun <- function(par, data_list, debug=0, ...){mmpp_ll(par, data_list, debug=0, ...)}
  } else{
    obj_fun <- function(par, data_list, debug=0, ...){mmpp_ll(par, data_list, debug=0, ...) - 2*pen_fun(par)}
  }
  
  if(debug==2) browser()
  
  # mmpp_ll(start, data_list, debug=1)
  
  if(fit){
    message('Optimizing likelihood...')  
    if(debug==2) browser()
    # opt <- nlminb(start=start, objective=mmpp_ll, data_list=data_list, ...)
    opt <- optimx::optimr(par=start, fn=obj_fun, method=method, data_list=data_list, ...)
    
    if(opt$convergence!=0){
      message("There was a problem with optimization... See output 'optimx' object.")
      # return(list(opt=opt, data_list=data_list))
      # hessian <- FALSE
      # V <- NULL
    }
    if(hessian){
      message('Calculating Hessian and variance-covariance matrices...')  
      H <- numDeriv::hessian(mmpp_ll, opt$par, data_list=data_list)
      V <- 2*solve(H)
    } else{
      V <- NULL
      H <- NULL 
    }
  } else{
    hessian <- FALSE
    V <- NULL
    opt <- list(par=start, objective=mmpp_ll(start, data_list))
  }
  
  if(debug==3) browser()
  
  ### Get real lambda values
  par <- as.vector(opt$par)
  # real <- get_reals(par, V, data_list, ddl, model_parameters)
  
  ### Get beta lambda values
  beta <- get_betas(par, V, data_list)
  reals <- get_reals(par, V, data_list, ddl, model_parameters)
  
  if(!hessian) V <- NULL
  
  out <- list(
    # par = c(beta_l,beta_q),
    par = par,
    vcov = V,
    log_lik = -0.5*opt$value,
    aic = opt$value + 2*length(par),
    results = list(
      beta = beta,
      real = reals
    ),
    opt = opt,
    start=start,
    data_list=data_list
  )
  
  if(debug==4) browser()
  
  return(out)
  
}