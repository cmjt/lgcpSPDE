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
struct covariate_list : vector<matrix <Type> > {
  covariate_list(SEXP x){  /* x = List passed from R */
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      (*this)(i) = asMatrix<Type>(sm);
    }
  }
};



template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density; // this where the structure for GMRF is defined
  using namespace Eigen;  // probably for sparseness calcs
  DATA_STRUCT(resp, counts_list);
  DATA_STRUCT(covariates, covariate_list); // list of design matrcies for each time step
  DATA_VECTOR(area);
  DATA_FACTOR(ID);
  DATA_STRUCT(spde,spde_t); // this as structure of spde object is defined in r-inla hence using that namespce
  PARAMETER_VECTOR(beta); // regression coefficients
  PARAMETER(log_kappa); // kappa of random field
  PARAMETER_ARRAY(x); //the random field/effect each matrix row is a time step
  PARAMETER(rho); // parameter of the AR(1) temporal process (only applicable for spatio-temporal model)
  int tsteps =  NLEVELS(ID); // number of time steps
  Type kappa = exp(log_kappa); // return the kappa parameter of the Random field
  Type nll = 0;
  SparseMatrix<Type> Q = Q_spde(spde,kappa); // create the precision matrix from the spde model for the GMRF
  //array<Type> st(resp(0).size(),tsteps);
  //st.setZero();
  nll = SEPARABLE(AR1(rho), GMRF(Q))(x); // x the random effect is a GMRF with precision Q
  for(int i = 0; i< (tsteps - 1); i++){
    vector<Type> eta = covariates(i)*beta; // design matrix and regression coeffs. (fixed effects)
    vector<Type> gmrf = (vector<Type> (x.col(i)));
    //std::cout << x.col(i).size() << "\n";
    vector<Type> respi = resp(i);
    for(int j = 0; j <respi.size(); j++){
      Type mu;
      mu = eta(j) + log(area(j)) + gmrf(j);
      Type lambda = exp(mu); // intensity 
      nll -= dpois(respi(j),lambda,true);
     
    }
  }
  ADREPORT(kappa);
  Type sigma2 = 1/(4*M_PI*pow(kappa,2));   // sigma2 = 1/(4*pi*k^2) with tau = 1
  ADREPORT(sigma2);
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  ADREPORT(range);
  return nll;
}


