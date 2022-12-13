library(moveMMPP)
library(Matrix)
library(tidyverse)
load("~/research/projects/methodology_devel/mmpp_movement/work/saved.RData")

# ddl$lambda$period_fac <- factor(ddl$lambda$period)

model_parameters <- list(
  lambda=list(form=~(island-1) + period:(island-1)),
  q=list(form=~from_mhi_zone_type-1)
)

cell_idx_df <- select(ddl$lambda, cell, cellx) %>% distinct()
data <- data %>% left_join(cell_idx_df, by="cell")

dml_list <- dm_lambda(model_parameters$lambda, ddl)
dmq_list <- dm_q(model_parameters$q, ddl)

data$period <- ifelse(is.na(data$cell), data$period-1, data$period)



N = as.integer(nrow(data))
ns = as.integer(length(unique(ddl$lambda$cell)))
np = as.integer(max(ddl$quad_pts$period))
# detection
id = data$idx-1
period = as.integer(data$period-1)
dt = data$delta
cell = as.integer(data$cellx-1)
# lambda
X_l = dml_list$X_l
off_l = dml_list$off_l
fix_l = dml_list$idx_l$fix
period_l = as.integer(dml_list$idx_l$period-1)
cell_l = as.integer(dml_list$idx_l$cell-1)
idx_l = as.integer(dml_list$idx_l$idx_l-1)
off_l = dml_list$off_l
# Q
from_q = as.integer(dmq_list$idx_q$from_cellx-1)
to_q = as.integer(dmq_list$idx_q$to_cellx-1)
X_q = dmq_list$X_q
off_q = dmq_list$off_q
idx_q = as.integer(dmq_list$idx_q$idx_q-1)

par_list <- list(
  beta_l=rep(0,ncol(dml_list$X_l)), 
  beta_q=rep(0, ncol(dmq_list$X_q))
)

beta_l <- par_list$beta_l
beta_q <- par_list$beta_q
Xb_l <- X_l %*% beta_l
Xb_q <- X_q %*% beta_q

from_to_q <- t(cbind(from_q, to_q))
Q <- moveMMPP:::load_Q(from_to_q, idx_q, Xb_q, off_q, ns)
L <- moveMMPP:::load_L(period_l, cell_l, idx_l, fix_l, Xb_l, off_l, ns, np)


st <- Sys.time()
moveMMPP:::mmpp_arma(id, period, dt, cell, ns, np, X_l, off_l, fix_l,
                     period_l, cell_l, idx_l, beta_l, from_to_q, X_q, 
                     off_q, idx_q, beta_q)
Sys.time()-st



QmLd <- ln_lik$QmLd
Q <- ln_lik$Q
L <- Matrix(diag(ln_lik$L_mat[,34]))

QmLd2 <- Q-L
max(abs(QmLd-QmLd2))



