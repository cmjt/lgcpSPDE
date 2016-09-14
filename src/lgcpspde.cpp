#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  DATA_VECTOR(points);
  DATA_VECTOR(area);
  DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hance using that namespce
  PARAMETER(beta0); // intercept term
  //PARAMETER(log_tau); 
  PARAMETER(log_kappa); // kappa of random field
  PARAMETER_VECTOR(x); //the random field/effect
  //Type tau = exp(log_tau); 
  Type kappa = exp(log_kappa); // return the kappa parameter of the Random field
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i <points.size(); i++){
    Type eta = beta0 + log(area(i)) + x(i); 
    Type lambda = exp(eta); // intensity 
    nll -= dpois(points(i),lambda,true); 
  }
  double nu = 1.0;            // nu = alpha-d/2 = 2-1 by eqn (2) in Lindgren 
  Type rho = sqrt(8*nu)/kappa;  // Distance at which correlation has dropped to 0.1 (p.  4 in Lindgren)
  ADREPORT(rho);
  ADREPORT(kappa);
  return nll;

}