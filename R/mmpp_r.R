#' @export
mmpp_r = function(par, data){
  N = data$N       # number of rows cell detection rows
  ns = data$ns       # number of sites
  np = data$np       # number of transition periods
  
  # DETECTIONS
  id = data$id        # individual id vector Nx1
  period = data$period+1   # vector of detection periods Nx1 for lambda_t model
  dt = data$dt        # time since last observation Nx1
  cell = data$cell+1     # observed cell detection Nx1
  cellNA = data$cellNA   # indicator of missing cell values.
  # LAMBDA
  X_l = data$X_l       # Design matrix for lambda
  fix_l = data$fix_l     # ns x np vector of fixed values for lambda
  period_l = data$period_l+1 # ns x np index for lambda periods (columns in Lmat)
  cell_l = data$cell_l+1  # ns x np index for lambda cells (rows in Lmat)
  idx_l = data$idx_l+1   # ns x np index for unique lambda values
  # Q MATRIX
  from_q = data$from_q+1   # index of 'from' cells for possible transitions
  to_q = data$to_q+1     # index of 'to' cells for possible transitions
  X_q = data$X_q      # Design matrix for Q
  idx_q = data$idx_q+1    # ns x ns index for unique Q values
  
  
  q_vals = exp(data$X_q %*% par$beta_q);
  l_vals = exp(data$X_l %*% par$beta_l);
  
  # Matrices for the HMM forward alg. 
  # matrix<Type> G(ns, ns);
  # matrix<Type> Q(ns,ns);
  # matrix<Type> QmLd(ns,ns);
  # matrix<Type> L_mat(ns,np);
  # vector<Type> L(ns);
  # vector<Type> P(ns);
  # matrix<Type> v(1,ns);
  # matrix<Type> phi(1,ns); phi.setZero();
  u = 0.0
  log_lik = 0
  log_lik_v = rep(0,N)
  
  
  Q = load_Q(from_q, to_q, idx_q, q_vals, ns)
  L_mat = load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np)
  # Start forward loop
  phi = matrix(0,1,ns)
  phi[1,cell[1]] = 1.0
  P = rep(0,ns)
  for(i in 2:N){
    # for(int i=1; i<1345; i++){
      if(id[i]!=id[i-1]){
        phi = 0.0*phi; phi[1,cell[i]] = 1
      } else{
        #Q = load_Q(from, to,Q_mat.col(period(i)), ns); #mod here for dynamic movement
        L = Matrix::Diagonal(x=L_mat[,period[i]])
        QmLd = as((Q-L)* dt[i], "dgCMatrix")
        #G = expm(QmLd);
        v = expQ2::v_exp_Q(phi, QmLd, 1.0e-10, renorm=FALSE, check=FALSE)
        # v =  expm::expAtv(t(QmLd), as.vector(phi))
        # v_m.row(i) = v.row(0);
        P = 0.0*P
        if(cellNA[i]==0){
          P[cell[i]] = L_mat[cell[i],period[i]];
          v = v * P #elementwise
        } 
        u = sum(v)
        log_lik_v[i] = log(u);
        phi = v/u;
      }
    } # end i
    return(-2*log_lik_v)
}