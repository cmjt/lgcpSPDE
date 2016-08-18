#' Function to fit a joint spatio-temporal model to geo-statistical data where one spatio-
#' temporal component is shared between the responses
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a matrix of response variables, each row corresponds to the spatial locations
#' in \code{locs}, the first column is the first response as a result of a
#' spatio-temporal process we want to ``copy'' into the second column response.
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T)
#' @param family a character vector specifying the assumed likelihoods of the response, by default is c("gaussian,"gaussian").
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{pcprior} is used with \code{param=c(0-0.9)}.
#' @param covariates a named data.frame of covariates
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#'
#' @export
geo.joint.fit <- function(mesh = NULL, locs.1 = NULL, locs.2 = NULL, response.1 = NULL, response.2 = NULL, t.index = NULL,  covariates = NULL, family = c("gaussian","gaussian"), verbose = FALSE, prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),hyper = list(theta=list(prior='normal', param=c(0,10))), control.inla=list(strategy='gaussian',int.strategy = 'eb')){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    # number of mesh nodes
    nv <- mesh$n
    if(!is.null(t.index)){
        temp <- t.index # temporal dimension
        k <- (mesh.t <- inla.mesh.1d(temp))$n # number of groups
        ## create projection matrix for loacations
        Ast.1 <- inla.spde.make.A(mesh = mesh, loc = locs.1, group = temp, n.group = k)
        Ast.2 <- inla.spde.make.A(mesh = mesh, loc = locs.2, group = temp, n.group = k)
        field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, group = temp, n.group = k)
        field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, group = temp, n.group = k)
        copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp, n.group = k)
        # temporal model "ar1"
        ctr.g <- list(model='ar1',param = prior.rho)
        if(!is.null(covariates)){
            m <- lgcpSPDE:::make.covs(covariates)
            cov.effetcs <- m[[1]]
            cov.form <- m[[2]]
            #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(response.1,NA)),
                                 A=list(Ast.1,1,1),
                                 effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)), cov.effets = cov.effects))
            formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                    " + f(field.1, model=spde, group = field.1.group, control.group=ctr.g)",
                    "+ f(field.2, model=spde, group = field.2.group , control.group=ctr.g)",
                    "+ f(copy.field, copy = field.1, fixed=FALSE )")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(response.1,NA)),
                                 A=list(Ast.1,1),
                                 effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
                  formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde, group = field.1.group, control.group=ctr.g) +
                      f(field.2, model=spde, group = field.2.group , control.group=ctr.g) +
                      f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
                 }
        stk.mark <- inla.stack(data=list(y=cbind(NA,response.2)),
                               A=list(Ast.2, Ast.2,1),
                               effects=list(field.2 = field.2, copy.field = copy.field, alpha0 = rep(1,nrow(locs.2))))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark)
    }else{
       
        ## create projection matrix for loacations
        Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1)
        Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2)
        field.1 <- field.2 <-  copy.field <-1:nv
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effects <- m[[1]]
            cov.form <- m[[2]]
                                        #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(response.1,NA)),
                                 A=list( Ast1,1,1),
                                 effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects))
            formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                    " + f(field.1, model=spde)",
                    "+ f(field.2, model=spde)",
                    "+ f(copy.field, copy = field.1, fixed=FALSE )")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(response.1,NA)),
                                 A=list( Ast1,1),
                                 effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
                  formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +
                      f(field.2, model=spde) +
                      f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
                 }
        stk.mark <- inla.stack(data=list(y=cbind(NA,response.2)),
                               A=list(Ast2, Ast2,1),
                               effects=list(field.2 = field.2, copy.field = copy.field, alpha0 = rep(1,nrow(locs.2))))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark)
        }
    ##call to inla
    result <- inla(as.formula(formula), family = family,
            data=inla.stack.data(stack),
            control.predictor=list(A=inla.stack.A(stack)),
            control.inla = control.inla,
            verbose = verbose,control.compute=list(dic=TRUE))
    result
}
                                      
    
