#' function to fit model used in bird paper
#' @return A \code{inla} result object
#'
#' @param locs BTO GBFS site locations
#' @param mesh Delauney triangulation of the UK
#' @param temp years of GBFS
#' @param binary.response list of length three. each referring to presence = 1, absence = 0 of the bird species
#' sparrowhawk, collared dove, house sparrow respectively
#' @param density.response list of length three. each referring to non-zero density of the bird species
#' sparrowhawk, collared dove, house sparrow respectively.
#' @param family a character vector of length two specifying the assumed likelihood of each species' response, by default
#' is rep(c("binomial","gamma"),3).
#' @param control.time (optional) supplied if the \code{temp} argument is given to fit a spatio-temporal model. This argument
#' controls the model and prior put on the hyperparameters of the model for the temporal component of the spatio-temporal
#' model. By default this is \code{list(model = 'ar1', param = list(theta = list(prior='pccor1', param = c(0, 0.9))))}
#' which is a pc.prior put on the rho coefficient of a AR(1) model with P(rho>0)=0.9.
#' @param control.inla a list which controls the fitting procedures INLA uses
#' by default this is \code{list(strategy='gaussian',int.strategy = 'eb')} for quick and dirty fitting.
#' @param hyper a list (of length 2) of lists of priors for each copy parameter. The first list has length 3
#' specifying the priors on the intra species interaction parameters ( c(beta_1, beta_2, beta_3) i.e., sparrowhawk, collared dove, house sparrow resp.).
#' the second element is of length 4 specifying the priors on the inter species interaction parameters
#' in order these refer to the parameters beta_z3, gamma_z3, beta_y3, and gamma_y3 (see model definition).
#' By default each is a N(0,10) (i.e.,  list(theta=list(prior='normal', param=c(0,10))))
#' @param control.compute a list of fit statistics the user wants INLA to return. By default this
#' is \code{list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE)}.
#' @param spde.new.params by default this is NULL. If supplied must be a named list with components:
#' \code{sig0} - typical standard deviation to use pc priors for hyperparams of spde model, \code{Psig} -
#' prob for sigma of pc prior, \code{rho0} - typical range to use pc priors for hyperparams of spde model,
#' and  \code{Prho} - prob for rho of pc prior
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
#' @param ... add inla options to speed up computation i.e., by giving starting values from a previous model
#'
#' @export

fit.multi <- function(locs = NULL, mesh = NULL, temp = NULL, binary.response = NULL, density.response = NULL,
                      family = rep(c("binomial","gamma"),3),
                      hyper = list(
                          intra =list(
                              beta_1 = list(theta=list(prior='normal', param=c(0,10))),
                              beta_2 = list(theta=list(prior='normal', param=c(0,10))),
                              beta_3 = list(theta=list(prior='normal', param=c(0,10)))),
                          inter = list(
                              beta_z3 = list(theta=list(prior='normal', param=c(0,10))),
                              gamma_z3 = list(theta=list(prior='normal', param=c(0,10))),
                              beta_y3 = list(theta=list(prior='normal', param=c(0,10))),
                              gamma_y3 = list(theta=list(prior='normal', param=c(0,10))))),
                      control.time = list(
                          model = 'ar1',param = list(
                                            theta = list(prior='pccor1', param = c(0, 0.9)))),
                      control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                      control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                      spde.new.params = NULL, verbose = FALSE, link = NULL,
                      ...){
    if(!is.null(spde.new.params)){
        rho0 <- spde.new.params$rho0
        Prho <- spde.new.params$Prho
        sig0 <- spde.new.params$sig0
        Psig <- spde.new.params$Psig
        spde <- lgcpSPDE:::inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    }else{
        spde <- inla.spde2.matern(mesh)
    }
    extra.args <- list(link = link)
    k <- length(table(temp))
    A <- inla.spde.make.A(mesh, loc = locs,group = temp)
    zsp <- binary.response[[1]]
    ysp <- density.response[[1]]
    zco <- binary.response[[2]]
    yco <- density.response[[2]]
    zhs <- binary.response[[3]]
    yhs <- density.response[[3]]
    rf1 <- inla.spde.make.index('rf1',n.spde=spde$n.spde,n.group=k,group = temp)
    brf1 <- inla.spde.make.index('brf1',n.spde=spde$n.spde,n.group=k,group = temp)
    rf2 <- inla.spde.make.index('rf2',n.spde=spde$n.spde,n.group=k,group = temp)
    brf2 <- inla.spde.make.index('brf2',n.spde=spde$n.spde,n.group=k,group = temp)
    g1rf1 <- inla.spde.make.index('g1rf1',n.spde=spde$n.spde,n.group=k,group = temp)
    g2rf2 <- inla.spde.make.index('g2rf2',n.spde=spde$n.spde,n.group=k,group = temp)
    rf3 <- inla.spde.make.index('rf3',n.spde=spde$n.spde,n.group=k,group = temp)
    g3rf1 <- inla.spde.make.index('g3rf1',n.spde=spde$n.spde,n.group=k,group = temp)
    g4rf2 <- inla.spde.make.index('g4rf2',n.spde=spde$n.spde,n.group=k,group = temp)
    brf3 <- inla.spde.make.index('brf3',n.spde=spde$n.spde,n.group=k,group = temp)
    b1 <- hyper$intra[[1]]
    b2 <- hyper$intra[[2]]
    b3 <- hyper$intra[[3]]
    g1 <- hyper[[2]]$inter[[1]]
    g2 <- hyper[[2]]$inter[[2]]
    g3 <- hyper[[2]]$inter[[3]]
    g4 <- hyper[[2]]$inter[[4]]
    stk.zsp <- inla.stack(data=list(y=cbind(zsp, NA,NA,NA,NA,NA)),
                          A=list(A), tag='sp.z',
                          effects=list(rf1 = rf1 ))
    stk.ysp<-inla.stack(data=list(y=cbind(NA,ysp,NA,NA,NA,NA)),tag="sp.y",
                        A=list(A),
                        effects=list(brf1 = brf1))
    stk.zco <- inla.stack(data=list(y=cbind(NA,NA,zco, NA,NA,NA)),
                          A=list(A), tag='co.z',
                          effects=list(rf2 = rf2))
    stk.yco<-inla.stack(data=list(y=cbind(NA,NA,NA,yco,NA,NA)),tag="co.y",
                        A=list(A),
                        effects=list(brf2 = brf2))
    stk.zhs <- inla.stack(data=list(y=cbind(NA,NA,NA,NA,zhs, NA)),
                          A=list(A,A,A,1), tag='hs.z',
                          effects=list(g1rf1 = g1rf1,g2rf2 = g2rf2,rf3 = rf3,
                                       alpha1=rep(1,length(zhs))))
    stk.yhs<-inla.stack(data=list(y=cbind(NA,NA,NA,NA,NA,yhs)),tag="hs.y",
                        A=list(A,A,A,1),
                        effects=list(g3rf1 = g3rf1,g4rf2 = g4rf2,brf3 = brf3,
                                     alpha2=rep(1,length(yhs))))
    stack<-inla.stack(stk.zsp,stk.ysp,stk.zco,stk.yco,stk.zhs,stk.yhs)
    formula<-y ~ 0 + alpha1 + alpha2 +
        f(rf1, model = spde,group = rf1.group, control.group = control.time) +
        f(brf1, copy = 'rf1' , fixed = FALSE, hyper = b1) +
        f(rf2, model=spde,group = rf2.group, control.group = control.time) +
        f(brf2, copy = 'rf2' , fixed = FALSE,hyper = b2) +
        f(g1rf1, copy = 'rf1' , fixed = FALSE,hyper = g1) +
        f(g2rf2, copy = 'rf2' , fixed = FALSE,hyper = g2) +
        f(g3rf1, copy = 'rf1' , fixed = FALSE,hyper = g3) +
        f(g4rf2, copy = 'rf2' , fixed = FALSE,hyper = g4) +
        f(rf3, model = spde,group = rf3.group,control.group = control.time) +
        f(brf3, copy = 'rf3' , fixed = FALSE,hyper = b3)
    result <- inla(as.formula(formula),
                   family = family,
                   data = inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack),link = extra.args$link),
                   control.inla = control.inla,
                   control.compute = control.compute,
                   verbose = verbose,
                   ...)
    
    
}
