#' Fitting a marked point process model
#'
#' Fits a spatio-temporal marked point process model
#' using INLA coupled with the SPDE approch 
#'
#' @return An inla model fit object, 
#' @param mesh delaunet triangulation of area, an object returned by \link{make.mesh} is suitable.
#' @param locs a matrix of \code{nrow} locations in \code{ncol} dimesions.
#' @param t.index a vector of length \code{nrow} of time index units refering to each point location
#' given in \link{locs}.
#' @param mark a vector of length \code{nrow} of marks refering to each point location
#' @param mark.family assumed likelihood for mark, by defalt "gaussian".
#' @param verbose Logical, if \code{TRUE}, model fitting is output
#' the console.
#'
#' @importMethodsFrom Matrix diag
#' @export
mark.pp.fit <- function(mesh = NULL, locs=NULL, t.index = NULL, mark = NULL, covariates = NULL, mark.family = "gaussian", verbose = FALSE, prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),hyper = list(theta=list(prior='normal', param=c(0,10)))){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    # number of observations
    n <- nrow(locs)
    # number of mesh nodes
    nv <- mesh$n
    if(!is.null(t.index)){
        temp <- t.index # temporal dimention
        k <- (mesh.t <- inla.mesh.1d(temp))$n # number of groups
        ## the response for the point pattern locations
        y.pp <- rep(0:1, c(k * nv, n))
        ## create projection matrix for loacations
        Ast <- inla.spde.make.A(mesh = mesh, loc = locs, group = temp, n.group = k)
        ##effect for LGCP used for point pattern
        st.volume <- diag(kronecker(Diagonal(n = k),spde$param.inla$M0))
        expected <- c(st.volume, rep(0, n))
        field.pp <- inla.spde.make.index('field.pp', n.spde = spde$n.spde, group = temp, n.group = k)
        field.mark <- inla.spde.make.index('field.mark', n.spde = spde$n.spde, group = temp, n.group = k)
        copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp, n.group = k)
        # temporal model "ar1"
        ctr.g <- list(model='ar1',param = prior.rho)
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effetcs <- m[[1]]
            cov.form <- m[[2]]
                                        #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA), e=expected),
                                 A=list(rBind(Diagonal(n=k*nv), Ast),1),
                                 effects=list(field.pp = field.pp, cov.effets = cov.effects))
            formula = paste("y", "~  0 + ", cov.form,
                    " + f(field.pp, model=spde, group = field.pp.group, control.group=ctr.g)",
                    "+ f(field.mark, model=spde, group = field.mark.group , control.group=ctr.g)",
                    "+ f(copy.field, copy = field.pp, fixed=FALSE )")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA), e=expected),
                                 A=list(rBind(Diagonal(n=k*nv), Ast)),
                                 effects=list(field.pp = field.pp))
                  formula <- y ~ 0 + f(field.pp, model=spde, group = field.pp.group, control.group=ctr.g) +
                      f(field.mark, model=spde, group = field.mark.group , control.group=ctr.g) +
                      f(copy.field, copy = "field.pp", fixed=FALSE, hyper = hyper )
                 }
        stk.mark <- inla.stack(data=list(y=cbind(NA,mark)),
                               A=list(Ast, Ast),
                               effects=list(field.mark = field.mark, copy.field = copy.field))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark)
    }else{
        y.pp <- rep(0:1, c( nv, n))
        ## create projection matrix for loacations
        Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
        ##effect for LGCP used for point pattern
        st.volume <- diag(spde$param.inla$M0)
        expected <- c(st.volume, rep(0, n))
                                        #fields
        field.pp <- field.mark <-  copy.field <-1:nv
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effetcs <- m[[1]]
            cov.form <- m[[2]]
                                        #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA), e=expected),
                                 A=list(rBind(Diagonal(n=nv), Ast),1),
                                 effects=list(field.pp = field.pp, cov.effets = cov.effects))
            formula = paste("y", "~  0 + ", cov.form,
                    " + f(field.pp, model=spde)",
                    "+ f(field.mark, model=spde)",
                    "+ f(copy.field, copy = field.pp, fixed=FALSE )")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA), e=expected),
                                 A=list(rBind(Diagonal(n=nv), Ast)),
                                 effects=list(field.pp = field.pp))
                  formula <- y ~ 0 + f(field.pp, model=spde) +
                      f(field.mark, model=spde) +
                      f(copy.field, copy = "field.pp", fixed=FALSE, hyper = hyper )
                 }
        stk.mark <- inla.stack(data=list(y=cbind(NA,mark)),
                               A=list(Ast, Ast),
                               effects=list(field.mark = field.mark, copy.field = copy.field))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark)
        }
    ##call to inla
    #family <-  dput(c("poisson",mark.family))
    result <- inla(formula, family = c("poisson",mark.family),
            data=inla.stack.data(stack),
            E=inla.stack.data(stack)$e,
            control.predictor=list(A=inla.stack.A(stack)),
            control.inla=list(strategy='gaussian',int.strategy = 'eb'),
            verbose = verbose)
    result
}
                                      
    
