
//** Test for NAs or non-finite values **//
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}
template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}

//** Structure for list of matrices from R **//
// template<class Type>
// struct MATRIX_LIST : vector< matrix<Type> >{
//   MATRIX_LIST(SEXP x){ // Constructor
//     (*this).resize(LENGTH(x));
//     for(int i=0; i<LENGTH(x); i++){
//       SEXP sm = VECTOR_ELT(x, i);
//       (*this)(i) = asMatrix<Type>(sm);
//     }
//   }
// };
// 
// //** Structure for list of vectors from R **//
// template<class Type>
// struct VECTOR_LIST : vector< vector<Type> >{
//   VECTOR_LIST(SEXP x){ // Constructor
//     (*this).resize(LENGTH(x));
//     for(int i=0; i<LENGTH(x); i++){
//       SEXP sm = VECTOR_ELT(x, i);
//       (*this)(i) = asVector<Type>(sm);
//     }
//   }
// };

//** Utility functions for HMM **//
template<class Type>
matrix<Type> mat_minus_diag(matrix<Type> M, vector<Type> d){
  matrix<Type> out = M;
  int n = M.rows();
  for(int i=0; i<n; i++){
    out(i,i) = M(i,i) - d(i);
  }
  return out;
}

template<class Type>
matrix<Type> load_Q(vector<int> from, vector<int> to, vector<int> idx_q, 
                    vector<Type> qvals, int ns){
  matrix<Type> Q(ns,ns); Q.setZero();
  int n = from.size();
  for(int i=0; i<n; i++){
    Q(from(i),to(i)) = qvals(idx_q(i));
  }
  return Q;
}

template<class Type>
matrix<Type> load_L(vector<int> period_l, vector<int> cell_l, 
                    vector<int> idx_l, vector<Type> fix_l, 
                    vector<Type> l_vals, int ns, int np){
  matrix<Type> L_mat(ns, np);
  int N = period_l.size();
  for(int i=0; i<N; i++){
    if(isNA(fix_l(i))){
      L_mat(cell_l(i), period_l(i)) = l_vals(idx_l(i));
    } else{
      L_mat(period_l(i), cell_l(i)) = fix_l(i);
    }
  }
  return L_mat;
}
