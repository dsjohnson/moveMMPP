#' @title Simulate a Markov Modulated Poisson Process
#' @param n Number of independent simulations
#' @param max_time ...
#' @param min_obs ...
#' @param Q A movement rate matrix
#' @param L A detection rate vector
#' @param start_loc Starting location of the path. If not specified, the limiting distribution 
#' implied by \code{Q} will be used to randomly select a starting location.
#' @importFrom Matrix diag
#' @importFrom stats rexp rpois runif
#' @author Devin S. Johnson
#' @export
sim_mmpp <- function(n=1, max_time, min_obs, Q, L, start_loc=NULL){
  out <- NULL
  for(i in 1:n){
    sim <- sim_mmpp_1(max_time, min_obs, Q, L, start_loc=NULL)
    sim$obs <- cbind(id=i, sim$obs)
    sim$path <- cbind(id=i, sim$path)
    out$obs <- rbind(out$obs, sim$obs)
    out$path <- rbind(out$path, sim$path)
  }
  return(out)
}

#' @importFrom Matrix diag
#' @importFrom stats rexp rpois runif
sim_mmpp_1 <- function(max_time, min_obs, Q, L, start_loc=NULL){
  # browser()
  ns <- length(L)
  exp_rate <- -Matrix::diag(Q)
  pi <- diag(1/exp_rate) %*% Q
  diag(pi) <- 0
  
  if(is.null(start_loc)){
    ud <- get_lim_ud(Q=Q)
    cur_loc <- sample.int(ns, 1, prob=ud)
  } else {
    cur_loc <- start_loc
  }
  cur_time <- 0
  path <- NULL
  path <- rbind(path, data.frame(timestamp=cur_time, loc=cur_loc))
  obs <- NULL
  total_obs <- 0
  stop_cond <- (cur_time >= max_time) | (total_obs >= min_obs)
  
  while(!stop_cond){
    move_time <- cur_time + rexp(1, rate=exp_rate[cur_loc])
    num_det <- rpois(1, L[cur_loc]*(move_time-cur_time))
    if(num_det > 0){
      det_times <- runif(num_det, cur_time, move_time)
      obs <- rbind(obs, data.frame(timestamp=sort(det_times), cell=cur_loc))
    }
    #
    cur_time <- move_time
    cur_loc <- sample.int(ns, 1, prob=pi[cur_loc,])
    path <- rbind(path, data.frame(timestamp=cur_time, loc=cur_loc))
    #
    total_obs <- total_obs + num_det
    stop_cond <- (cur_time >= max_time) | (total_obs >= min_obs)
  }
  st <- Sys.time()
  obs$timestamp <- st + (obs$timestamp)*86400
  path$timestamp <- st + (path$timestamp)*86400
  return(list(obs=obs, path=path))
  
  
}