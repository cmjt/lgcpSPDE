// modified the functions written by Skaug)
// 

namespace non_stat {
  using namespace Eigen;
  using namespace tmbutils;
  
  template<class Type>
  struct spde_nonstat_t{  
    SparseMatrix<Type> M0;        
    SparseMatrix<Type> M1;        
    SparseMatrix<Type> M2;
    matrix<Type> B0;
    matrix<Type> B1;
    matrix<Type> B2;
    spde_nonstat_t(SEXP x){  /* x = List passed from R */
      M0 = asSparseMatrix<Type>(getListElement(x,"M0"));
      M1 = asSparseMatrix<Type>(getListElement(x,"M1"));
      M2 = asSparseMatrix<Type>(getListElement(x,"M2"));
      B0 = asMatrix<Type>(getListElement(x,"B0"));
      B1 = asMatrix<Type>(getListElement(x,"B1"));
      B2 = asMatrix<Type>(getListElement(x,"B2"));
      
    }
  };
    
  template<class Type>
  SparseMatrix<Type> Q_spde_nonstat(spde_nonstat_t<Type> spde, Type theta1, Type theta2, Type theta3){
    vector<Type> log_phi0;
    vector<Type> log_phi1;
    vector<Type> log_phi2;
    vector<Type> log_phi12;
    log_phi0 = spde.B0.col(0) + spde.B0.col(1)*theta1 + spde.B0.col(2)*theta2 + spde.B0.col(3)*theta3;
    log_phi1 = spde.B1.col(0) + spde.B1.col(1)*theta1 + spde.B1.col(2)*theta2 + spde.B1.col(3)*theta3;
    log_phi2 = spde.B2.col(0) + spde.B2.col(1)*theta1 + spde.B2.col(2)*theta2 + spde.B2.col(3)*theta3;
    log_phi12 = log_phi1*log_phi2;
    SparseMatrix<Type> D0(log_phi0.size(),log_phi0.size());
    SparseMatrix<Type> D1(log_phi0.size(),log_phi0.size());
    SparseMatrix<Type> D2(log_phi0.size(),log_phi0.size());
    SparseMatrix<Type> D12(log_phi0.size(),log_phi0.size());
    std::cout << "length phi0 = " << log_phi1.size()  << "\n";
    for(int i=0; i <log_phi0.size(); i++){
      //std::cout << "i = " << i << "\n";
      for (int j=0; j<log_phi0.size(); j++){
    	if(i==j){ 
    	  D0.coeffRef(i,j) = exp(log_phi0(i));
    	  D1.coeffRef(i,j) = exp(log_phi1(i));
    	  D2.coeffRef(i,j) = exp(log_phi2(i));
    	  D12.coeffRef(i,j) = exp(log_phi12(i));
	  std::cout<< "& j = " << j << "\n";
    	}}}
    return D0*(D1*spde.M0*D1 + D12*spde.M1 + spde.M1.transpose()*D12 + spde.M2)*D0;
  }
} // end namespace non_stat
   
