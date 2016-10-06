#include <TMB.hpp>
#include "/home/charlotte/Documents/Gitwork/lgcpSPDE/inst/tmb/src/nonstatUN.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness cals
  using namespace non_stat_UN; // this where the structure for GMRF is defined my own)
  DATA_VECTOR(resp);
  DATA_VECTOR(area);
  DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_nonstat_UN_t); // this as structure defined by me
  PARAMETER(beta0); // intercept term
  PARAMETER_VECTOR(log_u); //diagonal of D1 matrix
  PARAMETER_VECTOR(x); //the random field/effect
  vector<Type> kappa = exp(log_u);
  SparseMatrix<Type> Q = Q_spde_nonstat_UN(spde,log_u); // create the precision matrix from the spde model for the GMRF
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i <resp.size(); i++){
    Type eta = beta0 + log(area(i)) + x(i); 
    Type lambda = exp(eta); // intensity 
    nll -= dpois(resp(i),lambda,true); // contribution from observed data
  }
  ADREPORT(kappa);
  return nll;

}
