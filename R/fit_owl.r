## function to fit owl marked pp model including snowfall information (that is covariate with misalignment)


fit.owl <- function(mesh = NULL, locs.pp=NULL, covariate = NULL, mark = NULL,
                        mark.family = "binomial",
                        verbose = FALSE,
                        hyper = list(theta=list(prior='normal', param=c(0,10))),
                        sig0 = 1,Psig = 0.5, rho0 = 0.3,Prho = 0.5,
                        control.inla = list(strategy='gaussian',int.strategy = 'eb'),
                        control.compute = list(dic = TRUE, waic = TRUE,cpo = TRUE, config = TRUE),
                        link = NULL, ...){
    spde <- lgcpSPDE:::inla.spde2.matern.new(mesh, prior.pc.rho = c(rho0, Prho), prior.pc.sig = c(sig0, Psig))
    extra.args <- list(link = link)
    m <- spde$n.spde
    n <- nrow(locs.pp)
    #projector matrix for locations of farms
    Ast <- inla.spde.make.A(mesh=mesh, loc=locs.pp)
    #pp location response
    y <- rep(0:1, c( m, n))
    #effort data for poisson pp
    st.volume <- diag(spde$param.inla$M0)
    expected <- c(st.volume, rep(0, n))
    A.pp<-rBind(Diagonal(n=m,rep(1,m)),Ast)
    
    stk.pp <- inla.stack(data=list(y=cbind(y,NA), e=expected),
                  A=list(A.pp),
                  tag="pp",
                  effects=list(list(field.pp=1:m)))
    #stack for response of interest
    stk.mark <- inla.stack(data=list(y=cbind(NA,mark)),
                  A=list( Ast,Ast,1),tag='mark',
                  effects=list(copy.field.pp = 1:m, mark.field = 1:m,
                               covariate = covariate))
    stack<-inla.stack(stk.pp,stk.mark)
    formula <- y ~ 0 + covariate + f(field.pp,model=spde) +
        f(copy.field.pp,copy="field.pp",hyper = hyper,fixed = FALSE) +
        + f(mark.field,model = spde)
    result <- inla(formula, family=c('poisson', mark.family),
            data=inla.stack.data(stack),
            control.predictor=list(A=inla.stack.A(stack),compute = TRUE,link = extra.args$link),
            control.inla = control.inla,
            control.compute = control.compute,
            verbose = verbose,
            ...)
}
