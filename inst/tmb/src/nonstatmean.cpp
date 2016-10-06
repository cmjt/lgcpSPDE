#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness clacs
  DATA_VECTOR(resp);
  DATA_VECTOR(area);
  DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespce
  DATA_MATRIX(dists); // Euclidean distance matrix of mesh nodes
  PARAMETER_VECTOR(beta0); // intercept term (spatially varying)
  PARAMETER(log_kappa); // kappa of random field
  PARAMETER(log_sigma); //parameter of the exponential function
  PARAMETER(rho); // parameter of the exponential covariance
  PARAMETER_VECTOR(x); //the random field/effect
  Type kappa = exp(log_kappa); // return the kappa parameter of the Random field
  Type sigma2 = exp(2.0*log_sigma); //return parameter of the exponential covariance for the expectation
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  matrix<Type> cov(resp.size(),resp.size()); //covariance matrix for spatially varying mean
  vector<Type> eta(resp.size()); //initialise intensity of Poisson
  eta = exp(log_sigma)*beta0;
  for(int i=0;i<resp.size();i++){
    cov(i,i)=Type(1); //ones on the diagonal
    for(int j=0;j<i;j++){
      cov(i,j)=exp(-rho*dists(i,j)); //exponential covariance of dists with parameter rho
      cov(j,i)=cov(i,j);
    }
  }
  MVNORM_t<Type> neg_log_density(cov); //define multivariate normal density with covariance cov  (exponential)
  Type nll = GMRF(Q)(x); // x the random effect is a GMRF with precision Q
  nll +=neg_log_density(beta0); //already defined as the negative log density
  for(int i = 0; i <resp.size(); i++){
    Type log_lam = eta(i) + log(area(i)) + x(i); 
    Type lambda = exp(log_lam); // intensity 
    nll -= dpois(resp(i),lambda,true); 
  }
  ADREPORT(kappa);
  ADREPORT(sigma2);
  return nll;
}
