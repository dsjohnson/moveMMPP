#' @title Fit MMPP movement model to resight data
#' @param data A processed data frame produced by the function \code{\link{process_data}}
#' @param ddl A design data list produced by the function \code{\link{make_design_data}}.
#' @param model_parameters Model formula for the detection and movement portions
#' of the MMPP model. Must be of the form \code{list(lambda=~form, q=~form)}.
#' @param hessian Logical. Should the Hessian matrix be calculated to obtain the parameter
#' variance-covariance matrix.
#' @param start Optional starting values for the parameter must be a list of the 
#' form \code{list(beta_l=c(), beta_q=c())}.
#' @param method Optimization method. See \code{\link[optimx]{optimr}}
#' @param ... Additional arguments passed to the optimization function 
#' \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#' @author Devin S. Johnson
#' @import optimx dplyr numDeriv
#' @export
fit_mmpp <- function(data, ddl, 
                     model_parameters=list(lambda=~1, q=~1), hessian=TRUE,
                     start=NULL, method="nlminb", ...
){
  
  cell <- cellx <- NULL
  
  cell_idx_df <- select(ddl$lambda, cell, cellx) %>% distinct()
  data <- data %>% left_join(cell_idx_df, by="cell")
  
  dml_list <- dm_lambda(model_parameters$lambda, ddl)
  dmq_list <- dm_q(model_parameters$q, ddl)
  
  data$period <- ifelse(is.na(data$cell), data$period-1, data$period)
  
  data_list <- list(
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
    par_list <- list(
      beta_l=rep(0,ncol(dml_list$X_l)), 
      beta_q=rep(0, ncol(dmq_list$X_q))
    )
  } else{
    par_list=start
  }
  
  start <- c(par_list$beta_l, par_list$beta_q)
  
  # message('Building model...')
  # foo <- MakeADFun(
  #   data=append(list(model="mmpp"), data_list),
  #   parameters=par_list,
  #   #random=c(),
  #   DLL="moveMMPP_TMBExports"
  # )
  
  message('Optimizing likelihood...')  
  # opt <- nlminb(start=start, objective=mmpp_ll, data_list=data_list, ...)
  opt <- optimx::optimr(par=start, fn=mmpp_ll, method=method, data_list=data_list, ...)
  
  if(opt$convergence!=0){
    message("There was a problem with optimization... See output 'optimx' object.")
    return(list(opt=opt, data_list=data_list))
  }
  if(hessian){
    message('Calculating Hessian and variance-covariance matrices...')  
    H <- numDeriv::hessian(mmpp_ll, opt$par, data_list=data_list)
    V <- 2*solve(H)
  }
  
  ### Get real lambda values
  X_l <- data_list$X_l
  lll <- 1:ncol(X_l)
  beta_l <- opt$par[lll]
  V_l <- V[lll,lll]
  l_vals <- exp(X_l %*% beta_l)
  # L <- load_L(data_list$period_l, data_list$cell_l, data_list$idx_l, 
  #                        data_list$fix_l, l_vals, data_list$ns, data_list$np)
  df_l <- bind_cols(
    select(ddl$lambda, cell, cellx, period, fix), 
    model.frame(model_parameters$lambda, ddl$lambda)
  )
  df_l$real <- l_vals[data_list$idx_l+1]
  df_l$real <- ifelse(is.na(df_l$real), df_l$fix, df_l$real)
  xbVbx_l <- X_l %*% V_l %*% t(X_l)
  se_real_l <- as.vector(l_vals) * sqrt(diag(xbVbx_l)) 
  ci_lower_l <- exp(X_l %*% beta_l - 1.96*sqrt(diag(V_l)))
  ci_upper_l <- exp(X_l %*% beta_l + 1.96*sqrt(diag(V_l)))
  df_l$se_real <- se_real_l[data_list$idx_l+1]
  df_l$se_real <- ifelse(is.na(df_l$se_real) & !is.na(df_l$fix), 0, df_l$se_real)
  df_l$ci_lower <- ci_lower_l[data_list$idx_l+1]
  df_l$ci_upper <- ci_upper_l[data_list$idx_l+1]
  df_l <- dplyr::distinct(df_l)
  
  ### Get beta lambda values
  l_nms <- colnames(X_l)
  df_beta_l <- data.frame(parameter = l_nms, est=beta_l, se_beta=diag(V_l))
  
  ### Get real Q  values
  X_q <- data_list$X_q
  qqq <- (ncol(X_l)+1):(ncol(X_l)+ncol(X_q))
  beta_q <- opt$par[qqq]
  V_q <- V[qqq,qqq]
  from_to_q <- t(cbind(data_list$from, data_list$to))
  q_vals <- exp(X_q %*% beta_q)
  q_nms <- colnames(X_q)
  # Q <- load_Q(from_to_q, data_list$idx_q, q_vals, data_list$ns)
  df_q <- model.frame(model_parameters$q, ddl$q)
  df_q$real <- q_vals[data_list$idx_q+1]
  xbVbx_q <- X_q %*% V_q %*% t(X_q)
  se_real_q <- as.vector(q_vals) * sqrt(diag(xbVbx_q)) 
  df_q$se_real <- se_real_q[data_list$idx_q+1]
  ci_lower_q <- exp(X_q %*% beta_q - 1.96*sqrt(diag(V_q)))
  ci_upper_q <- exp(X_q %*% beta_q + 1.96*sqrt(diag(V_q)))
  df_q$ci_lower <- ci_lower_q[data_list$idx_q+1]
  df_q$ci_upper <- ci_upper_q[data_list$idx_q+1]
  df_q <- dplyr::distinct(df_q)
  
  ### Get beta Q values
  q_nms <- colnames(X_q)
  df_beta_q <- data.frame(parameter = q_nms, est=beta_q, se_beta=diag(V_q))
  
  ####
  # statd <- eigen(t(Q))$vectors[,78] /sum(eigen(t(Q))$vectors[,78])
  # zones$ppp <- statd
  # mapview::mapview(zones, zcol='ppp')
  
  out <- list(
    par = c(beta_l,beta_q),
    vcov = V,
    log_lik = -0.5*opt$objective,
    results = list(
      beta = list(
        lambda = df_beta_l,
        q = df_beta_q
      ),
      real = list(
        lambda = df_l,
        q = df_q
      )
    ),
    opt = opt,
    start=start,
    data_list=data_list
  )
  
  return(out)
  
}