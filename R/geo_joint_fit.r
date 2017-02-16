#' Function to fit a  either a spatial or spatio-temporal joint model to geo-statistical data with an intercept and covariates for each likelihood
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a list of matrcies. The first element holds observation locations for the first likelihood, where each row corresponds to an observation. The second elemenr holds the observation locations for the second likelihood, each row corresponds to an observation. If no second element is supplied the observation locations for the first likelihood are used.
#' @param response a list (length two) of vectors of each response variable, each corresponds to the respective spatial locations
#' in \code{locs}.
#' @param temp (optional) a numeric vector specifying a temporal index for each observation (starting at 1.....T). It is assumed that the temporal indecies are common accross each reponse.
#' @param covariates (optional) a list (length 2) each element should contain a named data.frame of covariates. The first corresponding to the first likelihood, the second corresponding the the second likelihood.
#' @param family a character vector of length two specifying the assumed likelihood of each response, by default is c("gaussian","gaussian").
#' @param control.time (optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9. Assumed to be shared accross both responses
#' Refer to Simpson, martins, and rue for further details *****put in proper refs*****
#' @param control.inla a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param non.linear (optional) should be used if the user requires a non-linear covariate to be included in the model
#' Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
#' and \code{model} the random effect model the user wishes to use for \code{random.effect}
#' @param prediction (optional) should be used if the uses wnts to run a prediction model. Must be supplied as a
#' named list with \code{pred.locs} the locations where predictionas are to be made (referring to the first likelihood ONLY).
#' If a spatio-tempral model is to be fitted \code{pred.temp} the temporal indecies for the predictions (the same length as \code{pred.locs}).
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' 
#' @export

geo.joint.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian",
                    control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                    control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                    control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                    non.linear = NULL, prediction = NULL, verbose = FALSE){
     if(is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.j.fit(mesh = mesh, locs = locs, response = response,covariates = NULL,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.j.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(!is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.j.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp,
                                        family = family, covariates = NULL,
                                        control.time = control.time, control.inla = control.inla,
                                        control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(!is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.j.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                        temp = temp, family = family, control.time = control.time,
                                        control.inla = control.inla, control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    return(fit)
}





    
geo.spatial.j.temporal.fit <-function(mesh, locs, temp,...){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
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
            x = "\"field.1\""
            formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                    " + f(field.1, model=spde, group = field.1.group, control.group=ctr.g)",
                    "+ f(field.2, model=spde, group = field.2.group , control.group=ctr.g)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE )")
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
            x = "\"field.1\""
            formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                    " + f(field.1, model=spde)",
                    "+ f(field.2, model=spde)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE )")
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
                                      
    
