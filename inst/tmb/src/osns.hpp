// structure of oscillating precision of GRF with Q = k^4C +2k^2cos(pi theta)G GC^-1G
// 

namespace oscilate {
  using namespace Eigen;
  using namespace tmbutils;
  
  template<class Type>
  struct spde_oscilate_t{  
    SparseMatrix<Type> M0;      
    SparseMatrix<Type> M1;        
    SparseMatrix<Type> M2;
    spde_oscilate_t(SEXP x){  /* x = List passed from R */
      M0 = asSparseMatrix<Type>(getListElement(x,"M0"));
      M1 = asSparseMatrix<Type>(getListElement(x,"M1"));
      M2 = asSparseMatrix<Type>(getListElement(x,"M2"));
      
      
    }
  };
    
  template<class Type>
  SparseMatrix<Type> Q_spde_oscilate(spde_oscilate_t<Type> spde, Type kappa, Type theta){
    Type kappa_pow2 = kappa*kappa;
    Type kappa_pow4 = kappa_pow2*kappa_pow2;
    Type inlogittheta = exp(theta)/(1+exp(theta));       
    return kappa_pow4*spde.M0 + Type(2.0)*cos(M_PI*inlogittheta)*kappa_pow2*spde.M1 + spde.M2;  
  }
} // end namespace oscilate
   
