#' Function to fit a  spatio-temporal log-Gaussian Cox process model  with an intercept and covariates (optional)
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param mesh.pars a named vertor of mesh parameters, must contain
#' \code{cutoff} length at which to cut off triangle edge lengths,
#' \code{min} triangle edge length inside region,
#' and \code{max} triangle edge length outside region.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a vector of response variable, each corresponds to the spatial locations
#' in \code{locs}.
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T).
#' @param covariates a named data.frame of covariates 
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{pcprior} is used with \code{param=c(0-0.9)}. 
#' @param verbose Logical if \code{TRUE} model fit is output to screen.

#' @export


fit.lgcp <- function(mesh = NULL, mesh.pars = NULL, locs=NULL, temp = NULL, covariates = NULL, prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))), verbose = FALSE, control.inla=list(strategy='gaussian',int.strategy = 'eb'),return.attributes = FALSE){
    if(!is.null(covariates) & is.null(mesh)){
        stop("covariates must be supplied at the mesh nodes, thus, please supply mesh")
        }
    if(is.null(mesh)){
        if(is.null(mesh.pars)){
            warning("crude mesh constructed, highly recommended user supplies own triangulation")
            mesh <- make.mesh(locs = locs)
        }else{
            warning("highly recommended user checks triangulation")
            mesh <- make.mesh(loc = locs, mesh.pars = mesh.pars)
        }
    }else{
        mesh <- mesh
    }
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    # number of observations
    n <- nrow(locs)
    # number of mesh nodes
    nv <- mesh$n
    if(!is.null(temp)){
         k <- (mesh.t <- inla.mesh.1d(temp))$n
         Ast <- inla.spde.make.A(mesh=mesh, loc=locs, n.group=length(mesh.t$n),group=temp, group.mesh=mesh.t)
         field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
         st.volume <- diag(kronecker(Diagonal(n=k),spde$param.inla$M0))
         expected <- c(st.volume, rep(0, n))
         A.pp <- rBind(Diagonal(n=k *nv), Ast)
         ctr.g <- list(model='ar1',param = prior.rho)
         y.pp <- rep(0:1, c(k * nv, n))
         } else {
             Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
             field  <- 1:nv
             st.volume <- diag(spde$param.inla$M0)
             expected <- c(st.volume, rep(0, n))
             A.pp <- rBind(Diagonal(n=nv), Ast)
             y.pp <- rep(0:1, c(nv, n))
            }
    if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effects <- m[[1]]
            cov.form <- m[[2]]
            stack <- inla.stack(data=list(y=y.pp, e=expected),
                                 A=list(A.pp,1,1),
                                 effects=list(field = field, b0 = rep(1,length(y.pp)),cov.effects = cov.effects))
            if(!is.null(temp)){
                formula <- paste("y", "~  0  + b0 +", cov.form ,
                    " + f(field, model = spde, group = field.group, control.group = ctr.g)")
                }else{
                    formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
                }
            }else{
                if(!is.null(temp)){
                    stack <- inla.stack(data=list(y=y.pp, e=expected),
                                 A=list(A.pp,1),
                                 effects=list(field = field, b0 = rep(1,(k*nv)+n)))
                    formula <- y ~ 0 + b0 + f(field, model = spde, group = field.group, control.group = ctr.g)
                    }else{
                       stack <- inla.stack(data=list(y=y.pp, e=expected),
                                 A=list(A.pp,1),
                                 effects=list(field = field, b0 = rep(1,nv+n)))
                       formula <- y ~ 0  + b0 + f(field, model=spde)
                       }
                }
    ##call to inla
    result <- inla(as.formula(formula), family = "poisson",
            data=inla.stack.data(stack),
            E=inla.stack.data(stack)$e,
            control.predictor=list(A=inla.stack.A(stack)),
            control.inla = control.inla,
            verbose = verbose)
    if(return.attributes) attributes(result)$mesh <- as.list(mesh)
    result
}
                                      
