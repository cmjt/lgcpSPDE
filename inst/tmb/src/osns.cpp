#include <TMB.hpp>
#include "osns.hpp"


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  using namespace oscilate; // this where the structure for GMRF is defined my own)
  DATA_VECTOR(resp);
  DATA_VECTOR(area);
  DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_oscilate_t); // this as structure of spde object is defined in r-inla hance using that namespce
  PARAMETER(log_kappa);
  PARAMETER(theta); 
  PARAMETER(beta0); // intercept term
  PARAMETER_VECTOR(x); //the random field/effect
  Type kappa = exp(log_kappa);
  SparseMatrix<Type> Q = Q_spde_oscilate(spde,kappa,theta); // create the prec matrix from the spde model for the GMRF
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i <resp.size(); i++){
    Type eta = beta0 + log(area(i)) + x(i); 
    Type lambda = exp(eta); // intensity 
    nll -= dpois(resp(i),lambda,true); 
  }
  Type inlogittheta  = exp(theta)/(1+exp(theta)); 
  ADREPORT(inlogittheta);
  ADREPORT(kappa);
  return nll;
}
