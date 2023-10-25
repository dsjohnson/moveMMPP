Rcpp::sourceCpp("testing_info/arma_func_testing.cpp")


library(Matrix)
set.seed(98765)
n <- 10
# 5000 x 5000 matrices, 99% sparse
a <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)
b <- rsparsematrix(n, n, 0.01, rand.x=function(n) rpois(n, 1) + 1)
system.time(print(sum_by_iterator(a)))


library(Matrix)
set.seed(98765)
n <- 9
X <- rsparsematrix(n, n, 0.3)
diag(X) <- 1
X 
X2 <- sp_exp(X)
X22 <- sp_exp(X, row_sweep = FALSE)
apply(X2, 2, max)

X3 <- sp_exp(X)

max(X3)

X4 <- sp_norm(X3)


library(spatstat.sparse)
n <- 1000
X <- gridadjacencymatrix(c(n,n), across = TRUE, down = TRUE, diagonal=FALSE)
X <- as(X, "TsparseMatrix")
from_to <- t(cbind(X@i, X@j))
Xb_q_m <- rnorm(ncol(from_to))
Xb_q_r <- rep(0, nrow(X)) #rnorm(nrow(X))
eXb_q_r <- exp(Xb_q_r)
Q <- load_Q(from_to, Xb_q_r, Xb_q_m, ns=nrow(X))
# Q2 <- load_Q(from_to, Xb_q_r, Xb_q_m, ns=nrow(X), row_sweep = FALSE)


# library(rbenchmark)
# benchmark("V1" = {
#   load_Q(from_to, Xb_q_r, Xb_q_m, ns=nrow(X))
# },
# "V2" = {
#   load_Q(from_to, eXb_q_r, Xb_q_m, ns=nrow(X))
# },
# replications = 10,
# columns = c("test", "replications", "elapsed",
#             "relative", "user.self", "sys.self"))




