
// structure of the precision matrix of a non-stationary GRF with kappa=diag(Di) unknown
// 

namespace non_stat_UN {
  using namespace Eigen;
  using namespace tmbutils;
  
  template<class Type>
  struct spde_nonstat_UN_t{  
    SparseMatrix<Type> M0;      
    SparseMatrix<Type> M1;        
    SparseMatrix<Type> M2;
    matrix<Type> B2;
    spde_nonstat_UN_t(SEXP x){  /* x = List passed from R */
      M0 = asSparseMatrix<Type>(getListElement(x,"M0"));
      M1 = asSparseMatrix<Type>(getListElement(x,"M1"));
      M2 = asSparseMatrix<Type>(getListElement(x,"M2"));
      B2 = asMatrix<Type>(getListElement(x,"B2"));
      
    }
  };
    
  template<class Type>
  SparseMatrix<Type> Q_spde_nonstat_UN(spde_nonstat_UN_t<Type> spde, vector<Type> log_u){
    vector<Type> log_phi2;
    log_phi2 = spde.B2.col(0); // a column of 1s
    SparseMatrix<Type> D1(log_phi2.size(),log_phi2.size());
    SparseMatrix<Type> D2(log_phi2.size(),log_phi2.size());
    for(int i=0; i <log_phi2.size(); i++){
      for (int j=0; j<log_phi2.size(); j++){
	if(i==j){ 
	  D1.coeffRef(i,j) = exp(log_u(i));
	  D2.coeffRef(i,j) = exp(log_phi2(i));
	}}}
    return D1*spde.M0*D1 + D2*D1*spde.M1 + spde.M1.transpose()*D1*D2 + spde.M2;
  }
} // end namespace non_stat
   
