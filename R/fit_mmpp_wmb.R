#* @title Fit MMPP movement model to resight data using a batch approximation to ML
#* @param data A processed data frame produced by the function \code{\link{process_data}}
#* @param ddl A design data list produced by the function \code{\link{make_design_data}}.
#* @param model_parameters Model formula for the detection and movement portions
#* of the MMPP model. Must be of the form, e.g., 
#*  \code{list(lambda=list(form=~1, offset=NULL), q=list(form=~1, offset=~log(1.num_neigh)))}. For the 
#*  \code{q} model you must use \code{offset=0} to not have one. If it is left off, \code{offset=~log(1.num_neigh)))}
#*  will be used. 
#* @param batch_data A data set with columns labeled \code{id} and \code{batch} that describe which
#* which batch each individual belongs to for the weighted-mean batch ML approximation of 
#* Duncan (1980).
#* @param penalty_mat A covariance matrix for a multivariate normal penalty function. Using a 
#* penalized likelihood helps with optimization. If left as \code{penalty_mat=NULL} then
#* no penalty is placed on the likelihood. If each batch contains sufficient observations
#* to estimate all parameters, then no penalty is necessary. 
#* @param start Optional starting values for the parameter must be a list of the 
#* form \code{list(beta_l=c(), beta_q=c())}.
#* @param method Optimization method. See \code{\link[optimx]{optimr}}
#* @param fit Logical. Should the likelihood be optimized?
#* @param stop_opt_err Logical. Should the fitting stop (error) if any of the batches failes to converge
#* in the optimization. Otherwise a warning is thrown. 
#* @param verbose Messages/warnings for each batch fitting.
#* @param debug Integer from 1-4. Opens browser() at various points in the function call. Mostly for 
#* package developers. 
#* @param ... Additional arguments passed to the optimization function 
#* \code{\link[optimx]{optimr}} from the \code{\link[optimx]{optimx-package}}.
#* @author Devin S. Johnson
#* @references Duncan, G. M. (1980). Approximate maximum likelihood estimation 
#* with data sets that exceed computer limits. 
#* Journal of Econometrics, 14(2), 257-264.
#* @import optimx dplyr numDeriv foreach doFuture progressr
#* @importFrom stats ppois
#* @importFrom mvnfast dmvn
#* @export
#* 
# fit_mmpp_wmb <- function(data, ddl, 
#                          model_parameters = list(
#                            lambda = list(form=~1, offset=~NULL),
#                            q = list(form=~1, offset=~log(1/num_neigh)-1)
#                          ), 
#                          batch_data, penalty_mat=NULL, 
#                          start=NULL, method="nlminb", fit=TRUE,
#                          stop_opt_err = FALSE, verbose=TRUE, 
#                          debug=0, ...){
#   
#   cell <- cellx <- fix <- NULL
#   
#   if(debug==1) browser()
#   
#   cell_idx_df <- select(ddl$lambda, cell, cellx) %>% distinct()
#   data <- data %>% left_join(cell_idx_df, by="cell")
#   data <- data %>% left_join(batch_data, by='id')
#   
#   dml_list <- dm_lambda(model_parameters$lambda, ddl)
#   dmq_list <- dm_q(model_parameters$q, ddl)
#   
#   data$period <- ifelse(is.na(data$cell), data$period-1, data$period)
#   
#   data_list <- list(
#     N = as.integer(nrow(data)),
#     ns = as.integer(length(unique(ddl$lambda$cell))),
#     np = as.integer(max(ddl$quad_pts$period)),
#     # detection
#     id = data$idx-1,
#     batch = data$batch,
#     period = as.integer(data$period-1),
#     dt = data$delta,
#     cell = as.integer(data$cellx-1),
#     # lambda
#     X_l = dml_list$X_l,
#     off_l = dml_list$off_l,
#     fix_l = dml_list$idx_l$fix,
#     period_l = as.integer(dml_list$idx_l$period-1),
#     cell_l = as.integer(dml_list$idx_l$cell-1),
#     idx_l = as.integer(dml_list$idx_l$idx_l-1),
#     # Q
#     from_q = as.integer(dmq_list$idx_q$from_cellx-1),
#     to_q = as.integer(dmq_list$idx_q$to_cellx-1),
#     X_q = dmq_list$X_q,
#     off_q = dmq_list$off_q,
#     idx_q = as.integer(dmq_list$idx_q$idx_q-1)
#   )
#   
#   if(is.null(start)){
#     par_list <- list(
#       beta_l=rep(0,ncol(dml_list$X_l)), 
#       beta_q=rep(0, ncol(dmq_list$X_q))
#     )
#   } else{
#     par_list=start
#   }
#   
#   start <- c(par_list$beta_l, par_list$beta_q)
#   num_batch <- max(batch_data$batch)
#   if(!is.null(penalty_mat)){
#     cP <- chol(penalty_mat)
#     ln_pen <- function(par) (1/num_batch)*dmvn(par, rep(0,length(par)), cP, log=TRUE, isChol=TRUE)
#   } else {
#     ln_pen <- function(par) return(0)
#   }
#   
#   pen_n2ll <- function(par, data_list,...){mmpp_ll(par, data_list,...)-2*ln_pen(par)}
#   
#   opt_list <- vector("list", num_batch)
#   names(opt_list) <- paste0("batch_",1:num_batch)
#   
#   if(debug==2) browser()
#   
#   if(fit){
#     handlers(global=TRUE)
#     message("Beginning batch fits ...")
#     p <- progressor(num_batch)
#     opt_list <- foreach(b = c(1:num_batch), .options.future = list(seed = TRUE), .errorhandling="pass")%dofuture%{
#       batch_list <- batch_subset(data_list, b)
#       if(verbose) message(paste0('Optimizing likelihood for batch ', b, '...'))
#       # if(b>2) start <- rowMeans(sapply(opt_list[1:(b-1)], \(x)x$opt$par))
#       suppressWarnings(opt <- optimx::optimr(par=start, fn=pen_n2ll, method=method, data_list=batch_list, ...))
#       if(opt$convergence==0){
#         if(verbose) message(paste0('Calculating Hessian for batch ', b, '...'))
#         H <- numDeriv::hessian(pen_n2ll, opt$par, data_list=batch_list)/2
#         # V <- 2*solve(H)
#       } else{
#         H <- matrix(0, length(start), length(start))
#       }
#       p()
#       list(opt=opt, H=H)
#     }
#     cat("\n")
#     
#     if(any(sapply(opt_list, \(x) x$opt$convergence!=0))){
#       warning("There were fitting/optimization errors! Results may be unreliable. See returned 'optimr' list.") 
#     }
#   } else{
#     for(b in 1:num_batch){
#       opt <- list(par=start)
#       opt_list[[b]] <- list(opt=opt, H=diag(length(start)))
#     }
#   }
#   
#   if(debug==3) browser()
#   
#   message("Formatting output...")
#   
#   V <- Reduce("+", lapply(opt_list, \(x) x$H)) |> solve()
#   par <- as.vector(Reduce("+",lapply(opt_list, \(x) V%*%(x$H %*% x$opt$par))))
#   
#   if(!fit) V <- NULL
#   
#   real <- get_reals(par, V, data_list, ddl, model_parameters)
#   beta <- get_betas(par, V, data_list)
#   
#   opt <- list(
#     par = par,
#     value = mmpp_ll(par, data_list) -2*ln_pen(par),
#     opt_list=opt_list
#   )
#   
#   par_mat <- t(sapply(opt_list, \(x) x$opt$par))
#   colnames(par_mat) <- c(colnames(data_list$X_l), colnames(data_list$X_q))
#   
#   Vinv_list <- lapply(opt_list, \(x) x$H)
#   names(Vinv_list) <- paste0("batch_", c(1:length(Vinv_list)))
#   
#   out <- list(
#     par = par,
#     vcov = V,
#     log_lik = -0.5*opt$value,
#     results = list(
#       beta = beta,
#       real = real
#     ),
#     opt = opt,
#     start=start,
#     data_list=data_list,
#     penalty_mat = penalty_mat,
#     par_mat = t(sapply(opt_list, \(x) x$opt$par)),
#     Vinv_list = lapply(opt_list, \(x) x$H),
#     convergence = sapply(opt_list, \(x) x$opt$convergence==0)
#   )
#   
#   if(debug==4) browser()
#   
#   return(out)
# }
# 
# 
# batch_subset <- function(data_list, i){
#   data_list$cell <- data_list$cell[data_list$batch==i]
#   data_list$N <- length(data_list$cell)
#   data_list$id <- data_list$id[data_list$batch==i]
#   data_list$period <- data_list$period[data_list$batch==i]
#   data_list$dt <- data_list$dt[data_list$batch==i]
#   return(data_list)
# }