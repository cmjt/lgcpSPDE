#include <TMB.hpp>
// counts list (length one if model only spatial)
template<class Type>
struct counts_list : vector<vector <Type> >  {
  counts_list(SEXP x){  /* x = List passed from R for counts at each time step */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP counts = VECTOR_ELT(x, i);
      (*this)(i) = asVector<Type>(counts);
    }
  }
};


template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness calcs
  DATA_FACTOR(ID); // Factor ID of number of time steps (set to 1 for spatial only model)
  DATA_STRUCT(resp, counts_list);
  DATA_VECTOR(area);
  // DATA_IVECTOR(meshidxloc);
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespce
  PARAMETER(beta0); // intercept term
  PARAMETER(log_kappa); // kappa of random field
  PARAMETER_VECTOR(x); //the random field/effect
  PARAMETER(phi); // parameter of the AR(1) temporal process (only applicable for spatio-temporal model)
  int tsteps =  NLEVELS(ID); // number of time steps
  Type kappa = exp(log_kappa); // return the kappa parameter of the Random field
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  for(int j = 0; j < tsteps; j++){
    Type nll = SEPARABLE(AR1(phi)(x(1:respj.size() + (j - 1)*respj.size())),
			 GMRF(Q)(x(1:respj.size() + (j - 1)*respj.size()))); // x the random effect is a GMRF with precision Q
    vector<Type> respj = resp(j);
    for(int i = 0; i <respj.size(); i++){
      Type eta = beta0 + log(area(i)) + x(i); 
      Type lambda = exp(eta); // intensity 
      nll -= dpois(respj(i),lambda,true); 
    }
 } 
  ADREPORT(kappa);
  Type sigma2 = 1/(4*M_PI*pow(kappa,2));   // sigma2 = 1/(4*pi*k^2) with tau =1
  ADREPORT(sigma2);
  return nll;
}


  
