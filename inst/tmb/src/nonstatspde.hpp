// structure of non-stationary precision of GRF with kappa = diag(D1) = 
// exp(theta2 + thets3 (some user defined function of x and y)
// 

namespace non_stat {
  using namespace Eigen;
  using namespace tmbutils;
  
  template<class Type>
  struct spde_nonstat_t{  
    SparseMatrix<Type> M0;      
    SparseMatrix<Type> M1;        
    SparseMatrix<Type> M2;
    matrix<Type> B1;
    matrix<Type> B2;
    spde_nonstat_t(SEXP x){  /* x = List passed from R */
      M0 = asSparseMatrix<Type>(getListElement(x,"M0"));
      M1 = asSparseMatrix<Type>(getListElement(x,"M1"));
      M2 = asSparseMatrix<Type>(getListElement(x,"M2"));
      B1 = asMatrix<Type>(getListElement(x,"B1"));
      B2 = asMatrix<Type>(getListElement(x,"B2"));
      
    }
  };
    
  template<class Type>
  SparseMatrix<Type> Q_spde_nonstat(spde_nonstat_t<Type> spde, Type theta2, Type theta3){
    vector<Type> log_phi1;
    vector<Type> log_phi2;
    // the indecies should be 0,2,and 3 we are ignoring tau. Note ensure B1,B2 are nx4 matrices
    log_phi1 = spde.B1.col(0) +  spde.B1.col(2)*theta2 + spde.B1.col(3)*theta3;
    log_phi2 = spde.B2.col(0) +  spde.B2.col(2)*theta2 + spde.B2.col(3)*theta3;
    SparseMatrix<Type> D1(log_phi1.size(),log_phi1.size());
    SparseMatrix<Type> D2(log_phi1.size(),log_phi1.size());
    for(int i=0; i <log_phi1.size(); i++){
      for (int j=0; j<log_phi1.size(); j++){
	if(i==j){ 
	  D1.coeffRef(i,j) = exp(log_phi1(i));
	  D2.coeffRef(i,j) = exp(log_phi2(i));
	}}}
    return D1*spde.M0*D1 + D2*D1*spde.M1 + spde.M1.transpose()*D1*D2 + spde.M2;
  }
} // end namespace non_stat
   
