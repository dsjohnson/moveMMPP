#include <TMB.hpp>
// expm_generator:  v^T %*% exp(M) ->  M( from, to ) AND rowSums(M) = 1
template<class Type>
Eigen::SparseMatrix<Type> make_M( int CTMC_version,
                                  int n_g,
                                  matrix<int> At_zz,
                                  Type ln_D,
                                  vector<Type> h_g,
                                  vector<Type> colsumA_g ){

  int n_z = At_zz.rows();
  Type D = exp( ln_D );
  Eigen::SparseMatrix<Type> M_gg( n_g, n_g );
  if( CTMC_version==0 ){
    // Diffusion .. equal rate by cell, spread equally among neighbors
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += D / colsumA_g( At_zz(z,0) );
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= D / colsumA_g( At_zz(z,0) );
    }
    // Taxis
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += h_g(At_zz(z,1)) - h_g(At_zz(z,0));
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= h_g(At_zz(z,1)) - h_g(At_zz(z,0));
    }
  }else{
    // Combined taxis and diffusion
    for(int z=0; z<n_z; z++){
      M_gg.coeffRef( At_zz(z,0), At_zz(z,1) ) += D * exp( h_g(At_zz(z,1)) - h_g(At_zz(z,0)) );
      M_gg.coeffRef( At_zz(z,0), At_zz(z,0) ) -= D * exp( h_g(At_zz(z,1)) - h_g(At_zz(z,0)) );
    }
  }
  return M_gg;
}
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace sparse_matrix_exponential;

  // Data
  DATA_INTEGER( CTMC_version );
  DATA_INTEGER( Nmax );
  DATA_VECTOR( colsumA_g );
  DATA_MATRIX( L_gt );
  DATA_MATRIX( X_gz );
  DATA_IMATRIX( At_zz ); // Indices of A_gg
  DATA_INTEGER( gproj );
  //DATA_SPARSE_MATRIX( A_gg );
  //DATA_SPARSE_MATRIX( D_gg );

  // Parameters and random effects
  PARAMETER( ln_D );
  PARAMETER_VECTOR( gamma_z );

  // Global variables
  Type jnll = 0;
  int n_g = L_gt.rows();
  int n_t = L_gt.cols();
  vector<Type> nll_t( n_t );
  nll_t.setZero();
  matrix<Type> prob_gt( n_g, n_t );
  matrix<Type> pred_gt( n_g, n_t );
  matrix<Type> proj_gt( n_g, n_t );
  pred_gt.setZero();
  proj_gt.setZero();
  vector<Type> tmp_g( n_g );
  vector<Type> h_g(n_g);
  h_g = X_gz * gamma_z;
  Eigen::SparseMatrix<Type> Mrate_gg( n_g, n_g );
  //Mrate_gg = exp(ln_D)*A_gg - exp(ln_D)*D_gg;  // Diffusion-only case
  Mrate_gg = make_M( CTMC_version, n_g, At_zz, ln_D, h_g, colsumA_g );
  sparse_matrix_exponential::config<Type> myconfig = sparse_matrix_exponential::config<Type>();
  myconfig.Nmax = Nmax;
  expm_generator<Type> M_gg( Mrate_gg, myconfig );  // , sparse_matrix_exponential::config<Type>()
  
  // Project forward
  //Type renormalize = 0; 
  prob_gt.col(0) = L_gt.col(0);
  proj_gt(gproj,0) = Type(1.0); 
  for( int t=1; t<n_t; t++ ){
    // Predict movement (and re-normalize to correct for small numerical issues)
    pred_gt.col(t) = M_gg( vector<Type>(prob_gt.col(t-1)) ).array();
    pred_gt.col(t) = pred_gt.col(t) / pred_gt.col(t).sum();
    // Calculate probability (and re-normalize to avoid accumulating underflow in likelihood)
    prob_gt.col(t) = pred_gt.col(t).array() * L_gt.col(t).array();
    nll_t(t) = -1 * log( prob_gt.col(t).sum() );
    prob_gt.col(t) = prob_gt.col(t) / prob_gt.col(t).sum();
    // Project example release
    proj_gt.col(t) = M_gg( vector<Type>(proj_gt.col(t-1)) ).array();
  }
  jnll = nll_t.sum();
  
  // Reporting
  REPORT( Mrate_gg );
  REPORT( prob_gt );
  REPORT( pred_gt );
  REPORT( proj_gt );
  REPORT( nll_t );
  REPORT( h_g );
  return jnll;
}
