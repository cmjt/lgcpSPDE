
#' Function to fit a  either a spatial or spatio-temporal model to geo-statistical data with an intercept and covariates
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a vector of response variable, each corresponds to the spatial locations
#' in \code{locs}.
#' @param temp (optional) a numeric vector specifying a temporal index for each observation (starting at 1.....T).
#' @param covariates (optional) a named data.frame of covariates.
#' @param family a character vector specifying the assumed likelihood of the response, by default is "gaussian".
#' @param control.time (optional) supplied if the \code{temp} argumet is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9.
#' Refer to Simpson, martins, and rue for further details *****put in proper refs*****
#' @param control.inla a list which controls the fitting procedures INLA uses see Rue et al. ***ref book***
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param non.linear (optional) should be used if the user requires a non-linear covariate to be included in the model
#' Must be supplied as a named list with elements \code{random.effect} a numeric vector of the random effect indecies,
#' and \code{model} the random effect model the user wishes to use for \code{random.effect}
#' @param prediction (optional) should be used if the uses wnts to run a prediction model. Must be supplied as a
#' named list with \code{pred.locs} the locations where predictionas are to be made, and only if a spatio-tempral
#' model is to be fitted \code{pred.temp} the temporal indecies for the predictions (the same length as \code{pred.locs}).
#' @param sig0 by default = 1, typical standard deviation to use pc priors for hyperparams of spde model
#' @param rho0 by default = 0.3, typical range to use pc priors for hyperparams of spde model
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#'
#'
#'
#'
#' 
#' @export




geo.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian",
                    control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                    control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                    control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                    non.linear = NULL, prediction = NULL, sig0 = 1, rho0 = 0.3, verbose = FALSE){
    if(is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response,covariates = NULL,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction,sig0 = sig0, rho0 = rho0,
                               verbose = verbose)
    }
    if(is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, sig0 = sig0, rho0 = rho0,
                               verbose = verbose)
    }
    if(!is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp,
                                        family = family, covariates = NULL,
                                        control.time = control.time, control.inla = control.inla,
                                        control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction,sig0 = sig0, rho0 = rho0,
                                        verbose = verbose)
    }
    if(!is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                        temp = temp, family = family, control.time = control.time,
                                        control.inla = control.inla, control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, sig0 = sig0, rho0 = rho0,
                                        verbose = verbose)
    }
    return(fit)
}
    




#' spatial only model fitting
#'
geo.spatial.fit <- function(mesh, locs, response, covariates, family, control.inla, control.compute,
                            non.linear, prediction, sig0, rho0, verbose ){
    spde <-inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, 0.5), prior.pc.sig = c(sig0, 0.5))
    nv <- mesh$n
    n <- nrow(locs)
    if(is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n)))
        formula <- y ~ 0 + b0 +  f(field,model = spde)
    }
    if(is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        stk.obvs <- inla.stack(data=list(y = response),
                               A=list(index.o,1), tag = 'observation',
                               effects=list(field = 1:mesh$n, b0 = rep(1,n)))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0 +  f(field,model = spde)
    }
    if(!is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
    }
    if(!is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = response),
                            A=list(index.o, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field, model=spde)")
    }
    if(is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n),random = u))
        formula <- y ~ 0 + b0 +  f(random, model = u.mod) + f(field,model = spde)
    }
    if(!is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        index <- inla.spde.make.A(mesh = mesh, loc = locs)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = response),
                            A=list(index, 1, 1, 1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n),random = u, cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +  f(random, model = u.mod) +", cov.form," + f(field, model=spde)")
    }
    if(is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        stk.obvs <- inla.stack(data=list(y = response),
                               A=list(index.o,1,1), tag = 'observation',
                               effects=list(field = 1:mesh$n, b0 = rep(1,n), random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0  +  f(random, model = u.mod) +  f(field,model = spde)
    }
    if(!is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        locs.o <- locs
        locs.p <- prediction[["pred.locs"]]
        index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
        index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = response),
                            A=list(index.o, 1, 1,1), tag='observation',
                            effects=list(field = 1:mesh$n, b0 = rep(1,n), cov.effects = cov.effects, random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(index.p),
                              effects=list(field = 1:spde$n.spde), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0  +  f(random, model = u.mod) +", cov.form," + f(field, model=spde)")
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}






#' spatio temporal model fitting
#'
geo.spatial.temporal.fit <- function(mesh, locs, response, covariates, temp,  family,
                                    control.time, control.inla, control.compute,
                                    non.linear, prediction,sig0, rho0, verbose ){
    nv <- mesh$n
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    n <- nrow(locs)
    spde <-inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, 0.5), prior.pc.sig = c(sig0, 0.5))
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    field <- inla.spde.make.index('field',n.spde = spde$n.spde, group = temp, n.group = k)
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs ,group = temp, n.group = k)
    Y <- response
    if(is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n)))
        formula <- y ~ 0 + b0 +  f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        stk.obvs <- inla.stack(data=list(y = Y),
                               A=list(Ast,1), tag = 'observation',
                               effects=list(field = field, b0 = rep(1,n)))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0 +  f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&is.null(non.linear)&is.null(prediction)){
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +", cov.form,
                         " + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(!is.null(covariates)&is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0 +", cov.form," + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n),random = u))
        formula <- y ~ 0 + b0 +  f(random, model = u.mod) + f(field,model = spde,
                                                              group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&!is.null(non.linear)&is.null(prediction)){
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stack <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1, 1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n),random = u, cov.effects = cov.effects))
        formula <- paste("y", "~  0 + b0 +  f(random, model = u.mod) +", cov.form,
                         " + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    if(is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        stk.obvs <- inla.stack(data=list(y = Y),
                               A=list(Ast,1,1), tag = 'observation',
                               effects=list(field = field, b0 = rep(1,n), random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- y ~ 0 + b0  +  f(random, model = u.mod) +
            f(field,model = spde, group = field.group, control.group = control.time)
    }
    if(!is.null(covariates)&!is.null(non.linear)&!is.null(prediction)){
        locs.p <- prediction[["pred.locs"]]
        temp.p <- prediction[["pred.temp"]]
        k.p <- (mesh.t <- inla.mesh.1d(temp.p))$n
        Ast.p <- inla.spde.make.A(mesh = mesh, loc = locs.p ,group = temp.p, n.group = k.p)
        field.p <- inla.spde.make.index('field.p',n.spde = spde$n.spde, group = temp.p, n.group = k.p)
        u <- non.linear[["random.effect"]]
        u.mod <- non.linear[["model"]]
        m <- make.covs(covariates)
        cov.effects <- m[[1]]
        cov.form <- m[[2]]
        stk.obvs <- inla.stack(data=list(y = Y),
                            A=list(Ast, 1, 1,1), tag='observation',
                            effects=list(field = field, b0 = rep(1,n), cov.effects = cov.effects, random = u))
        stk.prd <- inla.stack(data=list(y = NA), A = list(Ast.p),
                              effects=list(field = field.p), tag='prediction')
        stack<-inla.stack(stk.obvs,stk.prd)
        formula <- paste("y", "~  0 + b0  +  f(random, model = u.mod) +",
                         cov.form," + f(field,model = spde, group = field.group, control.group = control.time)")
    }
    result <- inla(as.formula(formula), family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}



