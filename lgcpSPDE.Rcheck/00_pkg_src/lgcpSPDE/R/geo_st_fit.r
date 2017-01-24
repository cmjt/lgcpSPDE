#' Function to fit a  spatio-temporal model to geo-statistical data with an intercept and covariates
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a vector of response variable, each corresponds to the spatial locations
#' in \code{locs}.
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T).
#' @param covariates a named data.frame of covariates 
#' @param family a character vector specifying the assumed likelihood of the response, by default is "gaussian".
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{pcprior} is used with \code{param=c(0-0.9)}. 
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
geo.st.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian", prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))), verbose = FALSE){
    # number of mesh nodes
    nv <- mesh$n
    ## number of years and temporal mesh
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    # number of observations
    n <- nrow(locs)
    # the spde model for the latent field
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
    # the indecies for the latent spatio-temporal field
    field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
    # the projection matrix
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs ,group = temp, n.group = k)
    # define variables
    Y <- response
    ## covariates
    n.covs <- ncol(covariates)
    for(i in 1:n.covs){
        assign(colnames(covariates)[i],covariates[,i],envir = .GlobalEnv)
    }
    effects = list(field = field, sapply(colnames(covariates),get,simplify = FALSE),
                                   b0 = rep(1,length = n))
    # temporal model "ar1"
    ctr.g <- list(model='ar1',param = prior.rho)
    # model formula
    formula = paste("y", "~  0 + b0 +", paste(colnames(covariates),collapse = " + "),
                    " + f(field, model = spde, group = field.group, control.group = ctr.g)")
    stack <- inla.stack(data = list(y = Y),
                            A = list(Ast,1, 1),
                            tag = "response",
                            effects = effects)
    result <- inla(as.formula(formula),
                   family = family,
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.compute = list(config = TRUE),
                   verbose = verbose)
    attributes <- list()
    attributes$stack <- stack
    attributes(result) <- c(attributes(result), attributes)
    result

}
