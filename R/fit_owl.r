## function to fit owl marked pp model including snowfall information (that is covariate with misalignment)

fit.owl <- function(mesh = NULL, y = NULL, locs = NULL, covariates = NULL, misaligned.cov = NULL,
                    misaligned.locs = NULL,y.family = "binomial",mis.cov.family = "gamma",
                    verbose = FALSE, hyper = list(theta=list(prior='normal', param=c(0,10))),
                    sig0 = 1,Psig = 0.5, rho0 = 0.3,Prho = 0.5,
                    control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                    control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                    link = NULL, ...){
    if(is.null(misaligned.cov)){
        fit <- fit.owl.no.misaligned(mesh = mesh, y = y, locs = locs,covariates = covariates,
                                     y.family = y.family,
                                     verbose = verbose, sig0 = sig0, Psig = Psig,
                                     rho0 = rho0, Prho = Prho,
                                     control.inla = control.inla,
                                     control.compute = control.compute,
                                         link = link,
                                     ...)
    }
    if(!is.null(misaligned.cov)){
        fit <- fit.owl.misaligned(mesh = mesh, y = y, locs = locs,covariates = covariates,
                                  misaligned.cov = misaligned.cov,
                                  misaligned.locs = misaligned.locs,
                                  y.family = y.family,
                                  mis.cov.family = mis.cov.family,
                                  verbose = verbose, hyper = hyper,
                                  sig0 = sig0, Psig = Psig,
                                  rho0 = rho0, Prho = Prho,
                                  control.inla = control.inla,
                                      control.compute = control.compute,
                                  link = link,
                                  ...)
    }
    return(fit)
}

fit.owl.no.misaligned <- function(mesh, y, locs, covariates, 
                        y.family, verbose,sig0,Psig, rho0,Prho,
                        control.inla,control.compute,link, ...){
    spde <- lgcpSPDE:::inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    extra.args <- list(link = link)
    nv <- spde$n.spde
    #projector matrix for locations of farms
    Ast <- inla.spde.make.A(mesh=mesh, loc=locs)
    m <- make.covs(covariates)
    cov.effects <- m[[1]]
    cov.form <- m[[2]]
    stack <- inla.stack(data=list(y = y),
                  A=list(Ast,1),tag='mark',
                  effects=list(field = 1:nv,
                               covariate = cov.effects))
    formula <- paste("y", "~  0  +", cov.form," + f(field, model=spde)")
    result <- inla(as.formula(formula), family = y.family,
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack),compute = TRUE,link = extra.args$link),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    result
}

fit.owl.misaligned <- function(mesh, y, locs, covariates, 
                               misaligned.cov, misaligned.locs, y.family, mis.cov.family,
                               verbose, hyper, sig0,Psig, rho0,Prho,
                               control.inla,control.compute,link, ...){
    spde <- lgcpSPDE:::inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    extra.args <- list(link = link)
    nv <- spde$n.spde
    #projector matrix for locations of farms
    Ast <- inla.spde.make.A(mesh=mesh, loc=locs)
    Ast.mis <- inla.spde.make.A(mesh=mesh, loc = misaligned.locs)
    m <- make.covs(covariates)
    cov.effects <- m[[1]]
    cov.form <- m[[2]]
    stack.y <- inla.stack(data=list(y = cbind(y,NA)),
                  A=list(Ast,Ast,1),tag='mark',
                  effects=list(field = 1:nv,
                               copy.field = 1:nv,
                               covariate = cov.effects))
    stack.mis <- inla.stack(data=list(y = cbind(NA,misaligned.cov)),
                            A=list(Ast.mis),tag='misaligned',
                            effects=list(field.mis = 1:nv))
    stack <- inla.stack(stack.y, stack.mis)
     x = "\"field.mis\""
    formula <- paste("y", "~  0  +", cov.form," + f(field, model = spde) + f(field.mis, model = spde)",
                     "+ f(copy.field, copy =", x, ",fixed=FALSE, hyper = hyper )")
    result <- inla(as.formula(formula), family = c(y.family,mis.cov.family),
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack),compute = TRUE,link = extra.args$link),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    result
}
