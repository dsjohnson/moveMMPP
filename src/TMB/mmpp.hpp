/// @file mmpp.hpp

#include "include/helper.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// negative log-likelihood of the gamma distribution
template<class Type>
Type mmpp(objective_function<Type>* obj) {
  
  ////**** DATA ****////
  
  DATA_INTEGER(N);        // number of rows cell detection rows
  DATA_INTEGER(ns);       // number of sites
  DATA_INTEGER(np);       // number of transition periods
  
  // DETECTIONS
  DATA_VECTOR(id);        // individual id vector Nx1
  DATA_IVECTOR(period);   // vector of detection periods Nx1 for lambda_t model
  DATA_VECTOR(dt);        // time since last observation Nx1
  DATA_IVECTOR(cell);     // observed cell detection Nx1
  //LAMBDA
  DATA_MATRIX(X_l);       // Design matrix for lambda
  DATA_VECTOR(fix_l);     // ns x np vector of fixed values for lambda
  DATA_IVECTOR(period_l); // ns x np index for lambda periods (columns in Lmat)
  DATA_IVECTOR(cell_l);   // ns x np index for lambda cells (rows in Lmat)
  DATA_IVECTOR(idx_l);    // ns x np index for unique lambda values
  // Q MATRIX
  DATA_IVECTOR(from_q);   // index of 'from' cells for possible transitions
  DATA_IVECTOR(to_q);     // index of 'to' cells for possible transitions
  DATA_MATRIX(X_q);       // Design matrix for Q
  DATA_IVECTOR(idx_q);    // ns x ns index for unique Q values
  //DATA_VECTOR(p_q);     // vector of detection periods Nx1 for q_ij,t model
  //DATA_VECTOR(fix_q);   // vector of fixed values for Q;
  
  ////**** PARAMETERS ****////
  PARAMETER_VECTOR(beta_l); // Lambda parameters
  PARAMETER_VECTOR(beta_q); // movement parameters
  
  vector<Type> q_vals = exp(X_q * beta_q);
  vector<Type> l_vals = exp(X_q * beta_q);
  
  //int N = det.size();
  matrix<Type> G(ns, ns);
  matrix<Type> Q(ns,ns); 
  matrix<Type> QmLd(ns,ns);
  matrix<Type> L_mat(ns,np);
  vector<Type> L(ns);
  array<Type> P(ns);

  
  matrix<Type> v(1,ns);
  matrix<Type> phi(1,ns); phi.setZero();
  Type u = 0.0;
  Type log_lik = 0.0;
  
  Q = load_Q(from_q, to_q, idx_q, q_vals, ns);
  L_mat = load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np)
  
  // Start forward loop
  phi(0,cell(0)) = Type(1.0);
  for(int i=1; i<N; i++){
    
    if(id(i)!=id(i-1)){
      phi.setZero(); phi(0,cell(i)) = Type(1.0);
    } else{
      //Q = load_Q(from, to, qvals, ns); //mod here for dynamic movement
      L = L_mat.col(period(i));
      QmLd = mat_minus_diag(Q, L) * dt(i-1);
      G = expm(QmLd);
      v = phi * G;
      if(!isNA(cell(i))){
        P.setZero();
        P(cell(i)) = L(cell(i));
        v = v.array() * P; //elementwise
      }
      u = v.sum();
      log_lik += u;
      phi = v/u;
    }
    
  } // end i
  
  return -Type(2)*log_lik;
  
  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
