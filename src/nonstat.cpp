#include <TMB.hpp>
#include "nonstatspde.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  using namespace non_stat;
  DATA_VECTOR(resp);
  DATA_VECTOR(area);
  DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_nonstat_t); // this as structure of spde object is defined in r-inla hance using that namespce
  PARAMETER(theta1);
  PARAMETER(theta2);
  PARAMETER(theta3); 
  PARAMETER(beta0); // intercept term
  PARAMETER_VECTOR(x); //the random field/effect
  SparseMatrix<Type> Q = Q_spde_nonstat(spde,theta1,theta2,theta3); // create the prec matrix from the spde model for the GMRF  
  std::cout<< "Q rows = " << Q.rows() << "& Q cols = " << Q.cols() << "\n";
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i <resp.size(); i++){
    Type eta = beta0 + log(area(i)) + x(i); 
    Type lambda = exp(eta); // intensity 
    nll -= dpois(resp(i),lambda,true); 
  }
  return nll;
}
