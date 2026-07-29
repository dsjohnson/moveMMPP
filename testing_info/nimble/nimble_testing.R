library(nimble)
library(moveMMPP)

# 1. Point NIMBLE to installed headers and shared library
header_path <- system.file("include", "moveMMPP.h", package = "moveMMPP")
so_path     <- system.file("libs", paste0("moveMMPP", .Platform$dynlib.ext), package = "moveMMPP")

nimbleOptions(makeFlags = paste0(
  "-I\"", system.file("include", package = "RcppArmadillo"), "\" ",
  "-I\"", system.file("include", package = "expQ2"), "\" ",
  "-I\"", system.file("include", package = "Rcpp"), "\""
))

# 2. Register C function
c_mmpp_loglik_nimble <- nimbleExternalCall(
  prototype = function(
    id = double(1), n_id = integer(0),
    period = double(1), n_period = integer(0),
    dt = double(1), n_dt = integer(0),
    cell = double(1), n_cell = integer(0),
    ns = integer(0), np = integer(0),
    Xb_l = double(1), n_Xb_l = integer(0),
    fix_l = double(1), n_fix_l = integer(0),
    period_l = double(1), n_period_l = integer(0),
    cell_l = double(1), n_cell_l = integer(0),
    from_to = double(2), n_from_to_rows = integer(0), n_from_to_cols = integer(0),
    Xb_q_r = double(1), n_Xb_q_r = integer(0),
    Xb_q_m = double(1), n_Xb_q_m = integer(0),
    eq_prec = double(0),
    link_l = integer(0), link_r = integer(0), link_m = integer(0), struc = integer(0),
    a_l = double(0), a_r = double(0), a_m = double(0), norm_int = integer(0),
    log_p = integer(0)
  ) {},
  returnType = double(0),
  Cfun = "c_mmpp_loglik",
  headerFile = header_path,
  oFile = so_path
)

# 3. Define dMMPP density
dMMPP <- nimbleFunction(
  run = function(
    x = double(1), id = double(1), period = double(1), dt = double(1),
    ns = integer(0), np = integer(0),
    Xb_l = double(1), fix_l = double(1), period_l = double(1), cell_l = double(1),
    from_to = double(2), Xb_q_r = double(1), Xb_q_m = double(1),
    eq_prec = double(0, default = 1e-8),
    link_l = integer(0, default = 1), link_r = integer(0, default = 1),
    link_m = integer(0, default = 1), struc = integer(0, default = 1),
    a_l = double(0, default = 1.0), a_r = double(0, default = 1.0),
    a_m = double(0, default = 1.0), norm = integer(0, default = 1),
    log = integer(0, default = 0)
  ) {
    returnType(double(0))
    
    n_cell <- length(x)
    n_id <- length(id)
    n_period <- length(period)
    n_dt <- length(dt)
    
    n_Xb_l <- length(Xb_l)
    n_fix_l <- length(fix_l)
    n_period_l <- length(period_l)
    n_cell_l <- length(cell_l)
    
    n_ft_r <- dim(from_to)[1]
    n_ft_c <- dim(from_to)[2]
    
    n_Xb_q_r <- length(Xb_q_r)
    n_Xb_q_m <- length(Xb_q_m)
    
    return(c_mmpp_loglik_nimble(
      id, n_id, period, n_period, dt, n_dt, x, n_cell,
      ns, np,
      Xb_l, n_Xb_l, fix_l, n_fix_l, period_l, n_period_l, cell_l, n_cell_l,
      from_to, n_ft_r, n_ft_c,
      Xb_q_r, n_Xb_q_r, Xb_q_m, n_Xb_q_m,
      eq_prec, link_l, link_r, link_m, struc,
      a_l, a_r, a_m, norm, log
    ))
  }
)

# 4. Dummy Sampler & Distribution Registration
rMMPP <- nimbleFunction(
  run = function(
    n = integer(0), id = double(1), period = double(1), dt = double(1),
    ns = integer(0), np = integer(0), Xb_l = double(1), fix_l = double(1),
    period_l = double(1), cell_l = double(1), from_to = double(2),
    Xb_q_r = double(1), Xb_q_m = double(1), eq_prec = double(0, default = 1e-8),
    link_l = integer(0, default = 1), link_r = integer(0, default = 1),
    link_m = integer(0, default = 1), struc = integer(0, default = 1),
    a_l = double(0, default = 1.0), a_r = double(0, default = 1.0),
    a_m = double(0, default = 1.0), norm = integer(0, default = 1)
  ) {
    returnType(double(1))
    stop("rMMPP direct simulation is not implemented.")
  }
)

registerDistributions(list(
  dMMPP = list(
    BUGSdist = "dMMPP(id, period, dt, ns, np, Xb_l, fix_l, period_l, cell_l, from_to, Xb_q_r, Xb_q_m, eq_prec, link_l, link_r, link_m, struc, a_l, a_r, a_m, norm)",
    types = c(
      'value = double(1)', 'id = double(1)', 'period = double(1)', 'dt = double(1)',
      'ns = integer(0)', 'np = integer(0)',
      'Xb_l = double(1)', 'fix_l = double(1)', 'period_l = double(1)', 'cell_l = double(1)',
      'from_to = double(2)', 'Xb_q_r = double(1)', 'Xb_q_m = double(1)',
      'eq_prec = double(0)', 'link_l = integer(0)', 'link_r = integer(0)',
      'link_m = integer(0)', 'struc = integer(0)', 'a_l = double(0)',
      'a_r = double(0)', 'a_m = double(0)', 'norm = integer(0)'
    ),
    discrete = TRUE
  )
))