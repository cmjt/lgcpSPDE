#' Function to fit a  either a spatial or spatio-temporal joint model to geo-statistical data with an intercept and covariates for each likelihood
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a list of matrcies. The first element holds observation locations for the first likelihood, where each row corresponds to an observation. The second elemenr holds the observation locations for the second likelihood, each row corresponds to an observation. If no second element is supplied the observation locations for the first likelihood are used.
#' @param response a list (length two) of vectors of each response variable, each corresponds to the respective spatial locations
#' in \code{locs}.
#' @param temp (optional) a list of  numeric vectors specifying the temporal indcies for each response respectively.
#' @param covariates (optional) a list (length 2) each element should contain a named data.frame of covariates. The first corresponding to the first likelihood, the second corresponding the the second likelihood.
#' @param family a character vector of length two specifying the assumed likelihood of each response, by default is c("gaussian","gaussian").
#' @param control.time (optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9. Assumed to be shared accross both responses
#' Refer to Simpson, martins, and rue for further details *****put in proper refs*****
#' @param control.inla a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param hyper prior for the copy parameter by default is a N(0,10) i.e.,  list(theta=list(prior='normal', param=c(0,10)))
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param non.linear (optional) a list of named lists should be used if the user requires a non-linear covariate to be included for each likelihood. (i.e., non.linear = list(list(random.effect = idx.1, model = "iid"),list(random.effect = idx.2, model = "iid")) if the user wnats a iid effect for some idx.1 for the first likelihood and another for idx.2 for the second)
#' Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
#' and \code{model} the random effect model the user wishes to use for \code{random.effect}
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' 
#' @export

geo.joint.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL,
                          family = c("gaussian","gaussian"),
                          control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                          control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                          hyper = list(theta=list(prior='normal', param=c(0,10))),
                          control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                          non.linear = NULL,  verbose = FALSE){
    if(is.null(temp)){
        fit <- geo.spatial.j.fit(mesh = mesh, locs = locs, response = response,covariates = NULL,
                                 family = family, control.inla = control.inla, control.compute = control.compute,
                                 hyper = hyper,
                                 verbose = verbose)
    }
    if(!is.null(temp)&is.null(non.linear)){
        fit <- geo.spatial.j.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp,
                                          family = family, covariates = NULL,
                                          control.time = control.time, control.inla = control.inla,
                                          hyper = hyper, control.compute = control.compute,
                                          verbose = verbose)
    }
    if(!is.null(temp)&!is.null(non.linear)){
        fit <- geo.spatial.j.nl.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                             temp = temp, family = family, control.time = control.time,
                                             control.inla = control.inla, hyper = hyper, control.compute = control.compute,
                                             non.linear = non.linear,  verbose = verbose)
    }
    return(fit)
}


#' spatial only fitting
#' 
geo.spatial.j.fit <- function(mesh, locs, response, covariates, family,  control.inla,
                              hyper, control.compute, verbose){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    nv <- mesh$n
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2)
    field.1 <- field.2 <-  copy.field <-1:nv
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                             A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                             A=list( Ast2,Ast2,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field, alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form.1, cov.form.2,
                        " + f(field.1, model=spde)",
                        "+ f(field.2, model=spde)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                             A=list( Ast1,1),
                             effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                             A=list( Ast2,Ast2,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2))))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)

}


#' spatio-temporal model fitting 
geo.spatial.j.temporal.fit <-function(mesh, locs, response, covariates, temp, family, control.time, control.inla,
                              hyper, control.compute,  verbose){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    nv <- mesh$n
    temp <- temp # temporal dimension
    temp.1 <- temp[[1]]
    temp.2 <- temp[[2]]
    k.1 <- (mesh.t1 <- inla.mesh.1d(temp.1))$n
    k.2 <- (mesh.t2 <- inla.mesh.1d(temp.2))$n
    ## create projection matrix for loacations
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1, group = temp.1, n.group = k.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2, group = temp.2, n.group = k.2)
    field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, group = temp, n.group = k.1)
    field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, group = temp, n.group = k.2)
    copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp.2, n.group = k.2)
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field,
                                         alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        cov.form <- paste(cov.form.1,"+",cov.form.2)
        formula = paste("y", "~  0 + beta0 + alpha0 +", cov.form,
                        " + f(field.1, model=spde, group = field.1.group, control.group = control.time)",
                        "+ f(field.2, model=spde, group = field.2.group , control.group = control.time)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1))))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2))))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}



geo.spatial.j.nl.temporal.fit <-function(mesh, locs, response, covariates, temp, family, control.time, control.inla,
                              hyper, non.linear, control.compute,  verbose){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    nv <- mesh$n
    temp <- temp # temporal dimension
    temp.1 <- temp[[1]]
    temp.2 <- temp[[2]]
    k.1 <- (mesh.t1 <- inla.mesh.1d(temp.1))$n
    k.2 <- (mesh.t2 <- inla.mesh.1d(temp.2))$n
    ## create projection matrix for loacations
    response.1 <- response[[1]]
    response.2 <- response[[2]]
    locs.1 <- locs[[1]]
    locs.2 <- locs[[2]]
    u <- non.linear[[1]][["random.effect"]]
    u.mod <- non.linear[[1]][["model"]]
    u.2 <- non.linear[[2]][["random.effect"]]
    u.mod.2 <- non.linear[[2]][["model"]]
    Ast1 <- inla.spde.make.A(mesh = mesh, loc = locs.1, group = temp.1, n.group = k.1)
    Ast2 <- inla.spde.make.A(mesh = mesh, loc = locs.2, group = temp.2, n.group = k.2)
    field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, group = temp, n.group = k.1)
    field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, group = temp, n.group = k.2)
    copy.field <- inla.spde.make.index('copy.field', n.spde = spde$n.spde, group = temp.2, n.group = k.2)
    if(!is.null(covariates)){
        m.1 <- make.covs(covariates[[1]])
        m.2 <- make.covs(covariates[[2]])
        cov.effects.1 <- m.1[[1]]
        cov.form.1 <- m.1[[2]]
        cov.effects.2 <- m.2[[1]]
        cov.form.2 <- m.2[[2]]
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),cov.effects = cov.effects.1,
                                         random = u))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1,1,1),
                            effects=list(field.2 = field.2,copy.field = copy.field,
                                         alpha0 = rep(1,nrow(locs.2)),cov.effects = cov.effects.2, random.2 = u.2))
        stack <- inla.stack(stk.1,stk.2)
        x = "\"field.1\""
        cov.form <- paste(cov.form.1,"+",cov.form.2)
        formula = paste("y", "~  0 + beta0 + alpha0 + f(random, model = u.mod) + f(random.2, model = u.mod.2) +", cov.form,
                        " + f(field.1, model=spde, group = field.1.group, control.group = control.time)",
                        "+ f(field.2, model=spde, group = field.2.group , control.group = control.time)",
                    "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    }else{
        stk.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                            A=list( Ast1,1,1),
                            effects=list(field.1 = field.1, beta0 = rep(1,nrow(locs.1)),
                                         random = u))
        stk.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                            A=list( Ast2,Ast2,1),
                            effects=list(field.2 = field.2, copy.field = copy.field,alpha0 = rep(1,nrow(locs.2))))
        stack <- inla.stack(stk.1,stk.2)
        formula <- y ~ 0 + beta0 + alpha0 + f(field.1, model=spde) +  f(random, model = u.mod) +
            f(field.2, model=spde) +
            f(copy.field, copy = "field.1", fixed=FALSE, hyper = hyper )
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}
                         
