load("~/research/projects/methodology_devel/mmpp_movement/work/saved.RData")

from_to <- t(cbind(tmb_data$from, tmb_data$to))
qvals <- exp(tmb_data$X_q %*% tmb_par$beta_q)
idx_q <- as.integer(tmb_data$idx_q)
ns <- as.integer(tmb_data$ns)

Q <- moveMMPP:::load_Q(from_to, idx_q, qvals, ns)

beta_l <- rnorm(ncol(tmb_data$X_l), 0, 1)
period_l <- tmb_data$period_l
cell_l <- tmb_data$cell_l
idx_l <- tmb_data$idx_l
fix_l <- tmb_data$fix_l
l_vals <- exp(tmb_data$X_l %*% beta_l)
np <- tmb_data$np
ns <- tmb_data$ns

L <- moveMMPP:::load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np)




