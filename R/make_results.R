get_reals <- function(par, V=NULL, data_list, ddl, model_parameters){
  have_V <- !is.null(V)
  
  ### Lambda reals
  X_l <- data_list$X_l
  lll <- 1:ncol(X_l)
  beta_l <- par[lll]
  if(have_V) V_l <- V[lll,lll]
  l_vals <- exp(X_l %*% beta_l)
  vars_l <- unique(c(
    c('cell', 'cellx', 'period', 'fix'), 
    all.vars(model_parameters$lambda$form)
  ))
  df_l <- ddl$lambda[,vars_l]
  df_l$real <- l_vals[data_list$idx_l+1]
  df_l$real <- ifelse(is.na(df_l$real), df_l$fix, df_l$real)
  
  if(have_V){
    xbVbx_l <- X_l %*% V_l %*% t(X_l)
    se_real_l <- as.vector(l_vals) * sqrt(diag(xbVbx_l)) 
    ci_lower_l <- exp(X_l %*% beta_l - 1.96*se_real_l)
    ci_upper_l <- exp(X_l %*% beta_l + 1.96*se_real_l)
    df_l$se_real <- se_real_l[data_list$idx_l+1]
    df_l$se_real <- ifelse(is.na(df_l$se_real) & !is.na(df_l$fix), 0, df_l$se_real)
    df_l$ci_lower <- ci_lower_l[data_list$idx_l+1]
    df_l$ci_upper <- ci_upper_l[data_list$idx_l+1]
  }
  df_l$prob_det <- ppois(0, df_l$real, lower.tail=FALSE)
  if(have_V){
    df_l$ci_det_prob_lower <- ppois(0, df_l$ci_lower, lower.tail=FALSE)
    df_l$ci_det_prob_upper <- ppois(0, df_l$ci_upper, lower.tail=FALSE)
  }
  df_l <- dplyr::distinct(df_l)
  
  ### Q reals
  X_q <- data_list$X_q
  qqq <- (ncol(X_l)+1):(ncol(X_l)+ncol(X_q))
  beta_q <- par[qqq]
  if(have_V) V_q <- V[qqq,qqq]
  # from_to_q <- t(cbind(data_list$from, data_list$to))
  q_vals <- exp(X_q %*% beta_q)
  # q_nms <- colnames(X_q)

  vars_q <- all.vars(model_parameters$q$form)
  df_q <- ddl$q[,vars_q]
  df_q$real <- q_vals[data_list$idx_q+1]
  if(have_V){
    xbVbx_q <- X_q %*% V_q %*% t(X_q)
    se_real_q <- as.vector(q_vals) * sqrt(diag(xbVbx_q)) 
    df_q$se_real <- se_real_q[data_list$idx_q+1]
    ci_lower_q <- exp(X_q %*% beta_q - 1.96*se_real_q)
    ci_upper_q <- exp(X_q %*% beta_q + 1.96*se_real_q)
    df_q$ci_lower <- ci_lower_q[data_list$idx_q+1]
    df_q$ci_upper <- ci_upper_q[data_list$idx_q+1]
  }
  df_q <- dplyr::distinct(df_q)
  
  return(list(lambda=df_l, q=df_q))
  
}



get_betas <- function(par, V=NULL, data_list){
  have_V <- !is.null(V)
  
  ### Get lambda beta values
  lll <- 1:ncol(data_list$X_l)
  beta_l <- par[lll]
  if(have_V) V_l <- V[lll,lll]
  l_nms <- colnames(data_list$X_l)
  if(have_V){
    df_beta_l <- data.frame(parameter = l_nms, est=beta_l, se_beta=diag(V_l))
  } else{
    df_beta_l <- data.frame(parameter = l_nms, est=beta_l)
  }
  
  ### Get beta Q values
  qqq <- (ncol(data_list$X_l)+1):(ncol(data_list$X_l)+ncol(data_list$X_q))
  beta_q <- par[qqq]
  if(have_V) V_q <- V[qqq,qqq]
  q_nms <- colnames(data_list$X_q)
  if(have_V){
    df_beta_q <- data.frame(parameter = q_nms, est=beta_q, se_beta=diag(V_q))
  } else{
    df_beta_q <- data.frame(parameter = q_nms, est=beta_q)
  }
  
  return(list(lambda=df_beta_l, q=df_beta_q))
  
}