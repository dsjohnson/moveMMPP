library(moveMMPP)
library(Matrix)
library(tidyverse)
load("~/research/projects/methodology_devel/mmpp_movement/work/saved.RData")

ddl$lambda$period_fac <- factor(ddl$lambda$period)

model_parameters <- list(lambda=~(island - 1), 
                         q=~from_mhi_zone_type)

from_to_q <- t(cbind(tmb_data$from, tmb_data$to))
idx_q <- as.integer(tmb_data$idx_q)
beta_l <- rnorm(ncol(tmb_data$X_l), 0, 1)
beta_q <- rnorm(ncol(tmb_data$X_q), 0, 1)
period_l <- tmb_data$period_l
cell_l <- tmb_data$cell_l
idx_l <- tmb_data$idx_l
fix_l <- tmb_data$fix_l
X_l <- tmb_data$X_l
X_q <- tmb_data$X_q
id <- tmb_data$id
period <- tmb_data$period
dt <- tmb_data$dt
cell <- tmb_data$cell
np <- tmb_data$np
ns <- tmb_data$ns

l_vals <- exp(tmb_data$X_l %*% beta_l)
qvals <- exp(tmb_data$X_q %*% beta_q)
Q <- moveMMPP:::load_Q(from_to_q, idx_q, qvals, ns)
L <- moveMMPP:::load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np)


st <- Sys.time()
ln_lik <- moveMMPP:::mmpp_arma(id, period, dt, cell, ns, np, X_l, fix_l, period_l, 
          cell_l, idx_l, beta_l, from_to_q, X_q, idx_q, beta_q)
Sys.time()-st



QmLd <- ln_lik$QmLd
Q <- ln_lik$Q
L <- Matrix(diag(ln_lik$L_mat[,34]))

QmLd2 <- Q-L
max(abs(QmLd-QmLd2))



