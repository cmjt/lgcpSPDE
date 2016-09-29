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
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespce
  PARAMETER(beta0); // intercept term
  PARAMETER(log_kappa); // fixed kappa of random field
  PARAMETER_VECTOR(ns_log_kappa);
  PARAMETER(sigma_u); //variance on u (random effect on kappa
  vector<Type> u(points.size()); //the random bit of ns_log_kappa
  PARAMETER_VECTOR(x); //the random field/effect
  Type kappa = exp(log_kappa); // return the (fixed) kappa parameter of the Random field
  Type sum_u = 0;
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i <points.size(); i++){
    sum_u += u(i);
    ns_log_kappa(i) = log_kappa + sum_u;
    Type eta = beta0 + log(area(i)) + x(i); 
    Type lambda = exp(eta); // intensity 
    nll -= dpois(points(i),lambda,true); // contribution from observed data
    nll -=dnorm(u(i),Type(0),sigma_u,true);
  }
  ADREPORT(kappa);
  return nll;

}
