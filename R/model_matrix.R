#' @title Create sparse representations of design matrices for 
#' movement and resight models
#' @param par_list A list with named elements \code{form} ( the model formula) and 
#' \code{offset} (formula for the offset; must result in a single column from call to \code{model.matrix})
#' @param ddl Design data list
#' @details This function is not designed for end-users but is exported for posterity.
#' @name dm_matrix
NULL

#' @rdname dm_matrix
#' @importFrom stats model.frame model.matrix sd 
#' @export
dm_lambda <- function(par_list, ddl){
  X <- model.matrix(par_list$form, ddl$lambda)
  if(is.null(par_list$offset)){
    offset = rep(0,nrow(X))
  }else{
    offset = model.matrix(par_list$offset, ddl$lambda)
    if(ncol(offset)>1) stop("Lambda offset formula must result in model.matrix of only 1 column.")
    offset = as.vector(offset)
  }
  rX <- X[is.na(ddl$lambda$fix),]
  keep_col <- !((colMeans(rX)!=1) & (apply(rX, 2, sd)==0))
  X <- X[,keep_col]
  ### This section returns a reduce model matrix with only unique rows and 
  ### a lookup index. Right now this seems hard to rearrange and set up in 
  ### the TMB code. So, saving for later. 
  uX <- unique(X)
  uX <- uX[rowSums(uX)!=0,]
  dX <- data.frame(cbind(period=ddl$lambda$period, cellx=ddl$lambda$cellx, X))
  duX <- data.frame(cbind(idx_l=1:nrow(uX), uX))
  mX <- merge(dX, duX, all=TRUE)
  mX <- with(mX, mX[order(period, cellx),])
  lookup <- cbind(mX[,c('period','cellx','idx_l')], fix=ddl$lambda$fix)
  lookup$idx_l <- with(lookup, ifelse(!is.na(fix), NA, idx_l))
  return(list(X_l = uX, idx_l=lookup, off_l=offset))
}

#' @rdname dm_matrix
#' @export
dm_q <- function(par_list, ddl){
  X <- model.matrix(par_list$form, ddl$q)
  if(is.null(par_list$offset)){
    par_list$offset <- ~log(1/num_neigh)-1
    offset = model.matrix(par_list$offset, ddl$q)
    offset = as.vector(offset)
  } else if(par_list$offset==0){
    offset=rep(0,nrow(X))
  }else{
    offset = model.matrix(par_list$offset, ddl$q)
    if(ncol(offset)>1) stop("q offset formula must result in model.matrix of only 1 column.")
    offset = as.vector(offset)
    offset = ifelse(is.na(offset), 0, offset)
  }
  keep_col <- !((colMeans(X)!=1) & (apply(X, 2, sd)==0))
  X <- X[,keep_col]
  uX <- unique(X)
  uX <- uX[rowSums(uX)!=0,]
  dX <- data.frame(cbind(from_cellx=ddl$q$from_cellx, to_cellx=ddl$q$to_cellx, X))
  duX <- data.frame(cbind(idx_q=1:nrow(uX), uX))
  mX <- merge(dX, duX)
  mX <- with(mX, mX[order(from_cellx, to_cellx),])
  lookup <- mX[,c('from_cellx','to_cellx','idx_q')]
  return(list(X_q = uX, idx_q=lookup, off_q=offset))
}