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
  DATA_VECTOR(cellNA);   // indicator of missing cell values.
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
  vector<Type> l_vals = exp(X_l * beta_l);
  
  // Matrices for the HMM forward alg. 
  matrix<Type> G(ns, ns);
  matrix<Type> Q(ns,ns);
  matrix<Type> QmLd(ns,ns);
  matrix<Type> L_mat(ns,np);
  vector<Type> L(ns);
  //vector<Type> P(ns);
  matrix<Type> v(1,ns);
  matrix<Type> phi(1,ns); phi.setZero();
  Type u = 0.0;
  Type log_lik = Type(0);
  vector<Type> log_lik_v(N); log_lik_v.setZero();
  matrix<Type> phi_m(1345,ns);
  matrix<Type> v_m(1345,ns);
  
  
  Q = load_Q(from_q, to_q, idx_q, q_vals, ns);
  L_mat = load_L(period_l, cell_l, idx_l, fix_l, l_vals, ns, np);
  // 
  // // Start forward loop
  phi(0,cell(0)) = Type(1.0);
  // for(int i=1; i<N; i++){
  // for(int i=1; i<1345; i++){
  //   if(id(i)!=id(i-1)){
  //     phi.setZero(); phi(0,cell(i)) = Type(1.0);
  //   } else{
  //     //Q = load_Q(from, to,Q_mat.col(period(i)), ns); //mod here for dynamic movement
  //     L = L_mat.col(period(i));
  //     QmLd = mat_minus_diag(Q, L) * dt(i);
  //     G = expm(QmLd);
  //     v = phi * G;
  //     v_m.row(i) = v.row(0);
  //     if(cellNA(i)==0){
  //       //P.setZero();
  //       //P(cell(i)) = L(cell(i));
  //       //v = v.array() * P.array(); //elementwise
  //       //u = v.sum();
  //       log_lik += log(v(cell(i))*L(cell(i)));
  //       phi.setZero(); phi(0,cell(i)) = Type(1);
  //     } else{phi = v/v.sum();}
  //   }
  //   phi_m.row(i) = phi.row(0);
  // } // end i
  
  
  REPORT(Q);
  REPORT(phi_m);
  
  
  
  return -Type(2)*log_lik;
  
  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
