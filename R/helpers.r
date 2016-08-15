# takes in a named  matrix of covariates with ncol equal to the number of covariates
# returns a list containing the effects ready to be read by (inla.stack) and a covariate formula (ready to be read
# by inla
make.covs <- function(covariates){
    n.covs <- ncol(covariates)
    for(i in 1:n.covs){
        assign(colnames(covariates)[i],covariates[,i],envir = .GlobalEnv)
    }
    cov.effects <- sapply(colnames(covariates),get,simplify = FALSE)
    cov.form <- paste(colnames(covariates),collapse = " + ")
    return(list(cov.effects,cov.form))
}


## find which parts of the random field are inside the supplied spatial polygon
inwin<-function(proj, window){
    e<-expand.grid(proj$x,proj$y)
    o<-inside.owin(e[,1],e[,2],window)
    o<-matrix(o,nrow=length(proj$x))
}


## model fit for ant model (will generalise at some point) i.e. spatio temporal with two marks that is, p/a or ant nests, and p/a of pests

ants.pp.fit <- function(mesh = NULL, locs=NULL, marks = NULL, covariates = NULL, mark.family = c("binomial","binomial"), verbose = FALSE,hyper = list(theta=list(prior='normal', param=c(0,10)))){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    # number of observations
    n <- nrow(locs)
    # number of mesh nodes
    nv <- mesh$n
    ## the response for the point pattern locations
    y.pp <- rep(0:1, c( nv, n))
    mark.1 <- marks[,1]
    mark.2 <- marks[,2]
    ## create projection matrix for loacations
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
    ##effect for LGCP used for point pattern
    st.volume <- diag(spde$param.inla$M0)
    expected <- c(st.volume, rep(0, n))
    field.pp <- field.mark.1 <-  copy.field.mark.1 <- copy.field.1 <-  copy.field.2 <- 1:nv                            
    stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                         A=list(rBind(Diagonal(n=nv), Ast)),
                         effects=list(field.pp = field.pp))
    stk.mark.1 <- inla.stack(data=list(y=cbind(NA,mark.1,NA)),
                             A=list(Ast, Ast),
                             effects=list(field.mark.1 = field.mark.1, copy.field.1 = copy.field.1))
    stk.mark.2 <- inla.stack(data=list(y=cbind(NA,NA,mark.2)),
                             A=list(Ast, Ast ),
                             effects=list(copy.field.mark.1 = copy.field.mark.1, copy.field.2 = copy.field.2))
    stack <- inla.stack(stk.pp,stk.mark.1, stk.mark.2)
    formula <- y ~ 0 + f(field.pp, model=spde) +
        f(field.mark.1, model=spde) +
        f(copy.field.1, copy = "field.pp") +
        f(copy.field.2, copy = "field.pp" ) +
        f(copy.field.mark.1, copy = "field.mark.1", fixed=FALSE, hyper = hyper )
    ##call to inla
    result <- inla(formula, family = c("poisson",mark.family),
            data=inla.stack.data(stack),
            E=inla.stack.data(stack)$e,
            control.predictor=list(A=inla.stack.A(stack),compute = TRUE),
            control.inla=list(strategy='gaussian',int.strategy = "eb"),
            num.threads = 4,
            verbose = verbose)
    
    result
}
                                      
    

## function to fit owl marked pp model including snowfall information (that is covariate with misalignment) (need to generalise)
owls.pp.fit <- function(mesh = NULL, locs.pp=NULL, locs.cov = NULL, covariate = NULL,  mark = NULL,  mark.family = "binomial", covariate.family = "gaussian", verbose = FALSE, hyper = list(theta=list(prior='normal', param=c(0,10)))){
    spde <- inla.spde2.matern(mesh, alpha=2)
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
    #misalignet covariate data responses and prediction locs
    #covariate  projector matrix
    covAst<-inla.spde.make.A(mesh=mesh,loc=locs.cov)
    cov <- covariate
    #build data stacks
    stk.pp <- inla.stack(data=list(y=cbind(y,NA,NA), e=expected),
                  A=list(A.pp),
                  tag="pp",
                  effects=list(list(field.pp=1:m)))
    stk.cov <- inla.stack(data=list(y=cbind(NA,cov,NA)),
                  A=list(covAst),tag='covariate',
                  effects=list(field.cov = 1:m))
    #stack for response of interest
    stk.mark <- inla.stack(data=list(y=cbind(NA,NA,mark)),
                  A=list( Ast,Ast,Ast),tag='mark',
                  effects=list(copy.field.pp = 1:m, copy.field.cov = 1:m, mark.field = 1:m))
    stack<-inla.stack(stk.pp,stk.cov,stk.mark)
    formula <- y ~ 0  + f(field.pp,model=spde) + f(field.cov,model=spde) +
        f(copy.field.pp,copy="field.pp") +
        f(copy.field.cov,copy = "field.cov",fixed = FALSE, hyper = hyper) + f(mark.field,model = spde)
    #call to inla
    result <- inla(formula, family=c('poisson', covariate.family, mark.family),
            data=inla.stack.data(stack),
            control.predictor=list(A=inla.stack.A(stack),compute = TRUE),
            verbose = verbose)
    return(list(result=result, stack=stack))
}
    
                                       
    
