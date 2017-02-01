
#' Function to fit a  either a spatial or spatio-temporal model to geo-statistical data with an intercept and covariates
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
#'
#'
#'
#'
#' 
#' @export




geo.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian",
                    control.time = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                    control.inla=list(strategy='gaussian',int.strategy = 'eb'),
                    control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                    non.linear = NULL, prediction = NULL, verbose = FALSE){
    if(is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response, family = family,
                               control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                               family = family, control.inla = control.inla, control.compute = control.compute,
                               non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(!is.null(temp)&is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response, temp = temp, family = family,
                                        control.time = control.time, control.inla = control.inla,
                                        control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    if(!is.null(temp)&!is.null(covariates)){
        fit <- geo.spatial.temporal.fit(mesh = mesh, locs = locs, response = response,covariates = covariates,
                                        temp = temp, family = family, control.time = control.time,
                                        control.inla = control.inla, control.compute = control.compute,
                                        non.linear = non.linear, prediction = prediction, verbose = verbose)
    }
    return(fit)
}
    




#' spatial only model fitting
#'
geo.spatial.fit <- function(mesh, locs, response, covariates, family, control.inla, control.compute,
                            non.linear, prediction, verbose ){
    spde <- inla.spde2.matern(mesh, alpha = 2)
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
                   A=inla.stack.A(stack),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}






#' spatio temporal model fitting
#'
geo.spatial.temporalfit <- function(mesh, locs, response,temp, covariates, temp,  family,
                                    control.time, control.inla, control.compute,
                                    non.linear, prediction, verbose ){
    nv <- mesh$n
    k <- (mesh.t <- inla.mesh.1d(temp))$n
    n <- nrow(locs)
    spde <- inla.spde2.matern(mesh, alpha = 2)
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
            f(field,model = spde, group = field.group, control.group = control.time)y
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
                   A=inla.stack.A(stack),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose)
    return(result)
}
















geo.st.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL,covariates = NULL, family = "gaussian",
                       ctr.g = list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9)))),
                       control.inla=list(strategy='gaussian',int.strategy = 'eb'),
                       control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                       non-linear = NULL, prediction = NULL, verbose = FALSE){
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
                   control.inla  = control.inla,
                   verbose = verbose)
    attributes <- list()
    attributes$stack <- stack
    attributes(result) <- c(attributes(result), attributes)
    result

}





