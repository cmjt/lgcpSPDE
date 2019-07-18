#' Function to fit a  spatio-temporal log-Gaussian Cox process model  with an intercept and covariates (optional)
#'
#' @return A \code{INLA::inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{INLA::inla.mesh.2d}.
#' @param boundary spatial polygon of the point pattern observation window (optional). if supplied
#' weights at the mesh nodes outwith this will be set to zero.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T) (optional).
#' @param covariates a named data.frame of covariates (optional) 
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{INLA:::pcprior} is used with \code{param = c(0.9,0.9)}.
#' @param prior.range pc prior for the range of the latent field supplied as the vector c(range0,Prange)  (i.e., P(range < range0) = Prange), by default is \code{ c(5,0.9)}. NOTE should be changed to reflect range of the domain. 
#' @param prior.sigma pc prior for the sd of the latent field supplied as a vector (sigma0,Psigma) (i.e., P(sigma > sigma0) = Psigma), by default is \code{ c(1,0.005)}.
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param control.inla a list to control model fitting (as per inla)
#' @param control.fixed a list as per inla by default sets prior for precision intercept
#' @param ... other arguments taken by inla
#' @importMethodsFrom Matrix diag

#' @export


fit.lgcp <- function(mesh = NULL, boundary = NULL, locs = NULL, temp = NULL, covariates = NULL,
                     prior.rho = list(theta = list(prior='pccor1', param = c(0.0, 0.9))),
                     prior.range = c(5,0.9) ,
                     prior.sigma = c(1,0.005),
                     verbose = FALSE,
                     control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                     control.fixed = list(prec.intercept = 0.001),
                     return.attributes = FALSE,ns = NULL,
                     ...){
    mesh <- mesh
    if(!is.null(ns)){
        result <-
            fit.lgcp.ns(mesh = mesh,
                        locs = locs,
                        control.inla = control.inla,
                        verbose = verbose,
                        ns = ns)
    }else{
        spde <- inla.spde2.pcmatern(mesh = mesh,
                                    prior.range = prior.range,
                                    prior.sigma = prior.sigma)
        ## number of observations
        n <- nrow(locs)
        ## number of mesh nodes
        nv <- mesh$n
        if(!is.null(temp)){
            k <- (inla.mesh.1d(seq(1, k, by = 1)))$n
            Ast <- inla.spde.make.A(mesh = mesh,
                                    loc = locs,
                                    n.group = length(mesh.t$n),
                                    group = temp, group.mesh = mesh.t)
            field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
            if(!is.null(boundary)){
                w <- outwith(mesh = mesh, boundary = boundary)
                volume <- rep(w, k) * rep(diag(inla.mesh.fem(mesh.t)$c0), nv)
            }else{
                volume <- diag(kronecker(Diagonal(n = k),spde$param.inla$M0))
            }
            expected <- c(volume, rep(0, n))
            A.pp <- rBind(Diagonal(n = k *nv), Ast)
            ctr.g <- list(model = 'ar1',param = prior.rho)
            y.pp <- rep(0:1, c(k * nv, n))
        }else{
            Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
            field  <- 1:nv
            if(!is.null(boundary)){
                w <- outwith(mesh = mesh, boundary = boundary)
                volume <- w * diag(spde$param.inla$M0)
            }else{
                volume <- diag(spde$param.inla$M0)
            }
            expected <- c(volume, rep(0, n))
            A.pp <- rBind(Diagonal(n = nv), Ast)
            y.pp <- rep(0:1, c(nv, n))
        }
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effects <- m[[1]]
            cov.form <- m[[2]]
            stack <- inla.stack(data=list(y=y.pp, e=expected),
                                A=list(A.pp,1,1),
                                effects=list(field = field,
                                             b0 = rep(1,length(y.pp)),
                                             cov.effects = cov.effects))
            if(!is.null(temp)){
                formula <- paste("y", "~  0  + b0 +", cov.form ,
                                 " + f(field, model = spde, group = field.group,control.group = ctr.g)")
            }else{
                formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
            }
        }else{
            if(!is.null(temp)){
                stack <- inla.stack(data = list(y = y.pp, e = expected),
                                    A = list(A.pp,1),
                                    effects = list(field = field, b0 = rep(1,(k*nv)+n)))
                formula <- y ~ 0 + b0 + f(field, model = spde, group = field.group,
                                          control.group = ctr.g)
            }else{
                stack <- inla.stack(data = list(y = y.pp, e = expected),
                                    A = list(A.pp,1),
                                    effects = list(field = field, b0 = rep(1,nv+n)))
                formula <- y ~ 0  + b0 + f(field, model = spde)
            }
        }
        ##call to inla
        result <- inla(as.formula(formula), family = "poisson",
                       data = inla.stack.data(stack),
                       E = inla.stack.data(stack)$e,
                       control.predictor = list(A = inla.stack.A(stack),compute = TRUE),
                       control.inla = control.inla,
                       control.fixed = control.fixed,
                       verbose = verbose,
                       ...)
        if(return.attributes) {
            attributes(result)$mesh <- as.list(mesh)
            attributes(result)$volume <- volume
            }
    }
    result
}

