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


                                      
#############################
fit.two.marks.lgcp <- function(mesh = NULL, locs=NULL, t.index = NULL, mark.1 = NULL,  mark.2 = NULL, covariates = NULL, mark.family = c("gaussian","gaussian"), verbose = FALSE, prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),hyper = list(theta=list(prior='normal', param=c(0,10)))){
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
        cpoy.field.mark <- inla.spde.make.index('copy.field.mark', n.spde = spde$n.spde, group = temp, n.group = k)
        copy.field2 <- inla.spde.make.index('copy.field2', n.spde = spde$n.spde, group = temp, n.group = k)
        # temporal model "ar1"
        ctr.g <- list(model='ar1',param = prior.rho)
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effetcs <- m[[1]]
            cov.form <- m[[2]]
                                        #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                                 A=list(rBind(Diagonal(n=k*nv), Ast),1),
                                 effects=list(field.pp = field.pp, cov.effets = cov.effects))
            x = "\"field.pp\""
            y =  "\"field.mark\""
            formula = paste("y", "~  0 ", cov.form,
                    " + f(field.pp, model=spde, group = field.pp.group, control.group=ctr.g)",
                    "+ f(field.mark, model=spde, group = field.mark.group , control.group=ctr.g)",
                    "+ f(copy.field, copy =",x," )",
                    "+ f(copy.field2, copy =",x,")",
                    " + f(copy.field.mark, copy =",y,",fixed = FALSE, hyper = hyper)")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                                 A=list(rBind(Diagonal(n=k*nv), Ast)),
                                 effects=list(field.pp = field.pp))
                  formula <- y ~ 0 + f(field.pp, model=spde, group = field.pp.group, control.group=ctr.g) +
                      f(field.mark, model=spde, group = field.mark.group , control.group=ctr.g) +
                      f(copy.field, copy = "field.pp") +
                      f(copy.field2, copy = "field.pp") +
                      f(copy.field.mark, copy = "field.mark", fixed = FALSE, hyper = hyper )
                 }
        stk.mark.1 <- inla.stack(data=list(y=cbind(NA,mark.1,NA)),
                               A=list(Ast, Ast),
                               effects=list(field.mark = field.mark, copy.field = copy.field))
        stk.mark.2 <- inla.stack(data=list(y=cbind(NA,NA,mark.2)),
                               A=list(Ast, Ast),
                               effects=list(copy.field.mark = copy.field.mark, copy.field2 = copy.field2))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark.1,stk.mark.1)
    }else{
        y.pp <- rep(0:1, c( nv, n))
        ## create projection matrix for loacations
        Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
        ##effect for LGCP used for point pattern
        st.volume <- diag(spde$param.inla$M0)
        expected <- c(st.volume, rep(0, n))
                                        #fields
        field.pp <- field.mark <-  copy.field <- copy.field2 <- copy.field.mark <-1:nv
        if(!is.null(covariates)){
            m <- make.covs(covariates)
            cov.effetcs <- m[[1]]
            cov.form <- m[[2]]
                                        #create data stacks
            stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                                 A=list(rBind(Diagonal(n=nv), Ast),1),
                                 effects=list(field.pp = field.pp, cov.effets = cov.effects))
            x = "\"field.pp\""
            y =  "\"field.mark\""
            formula = paste("y", "~  0  + ", cov.form,
                    " + f(field.pp, model=spde)",
                    "+ f(field.mark, model=spde)",
                    "+ f(copy.field, copy =",x,")",
                    "+ f(copy.field2, copy =",x,")",
                    " + f(copy.field.mark, copy =",y,",fixed = FALSE, hyper = hyper)")
            }else{
                 stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                                 A=list(rBind(Diagonal(n=nv), Ast)),
                                 effects=list(field.pp = field.pp))
                  formula <- y ~ 0  + f(field.pp, model=spde) +
                      f(field.mark, model=spde) +
                      f(copy.field, copy = "field.pp")+
                      f(copy.field2, copy = "field.pp") +
                      f(copy.field.mark, copy = "field.mark", fixed = FALSE, hyper = hyper)
                 }
        stk.mark.1 <- inla.stack(data=list(y=cbind(NA,mark.1,NA)),
                               A=list(Ast, Ast),
                               effects=list(field.mark = field.mark, copy.field = copy.field))
        stk.mark.2 <- inla.stack(data=list(y=cbind(NA,NA,mark.2)),
                               A=list(Ast, Ast),
                               effects=list(copy.field.mark = copy.field.mark, copy.field2 = copy.field2))
        ## combine data stacks
        stack <- inla.stack(stk.pp,stk.mark.1,stk.mark.2)
        }
    ##call to inla
    result <- inla(as.formula(formula), family = c("poisson",mark.family),
            data=inla.stack.data(stack),
            E=inla.stack.data(stack)$e,
            control.predictor=list(A=inla.stack.A(stack)),
            control.inla=list(strategy='gaussian',int.strategy = 'eb'),
            verbose = verbose)
    result
}
                                      
    

#############################

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
    
                                       
######################### fit non standard i.e. non-stationary or TMB    
fit.ns.kappa.inla <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["fn"]]))stop("conjectured non-stationary function must be supplied")
    fn <- ns[["fn"]]
    B.kappa <- cbind(0,0,1,fn)
    spde <- inla.spde2.matern(mesh = mesh,
                              alpha = 2, B.tau = cbind(0,1,0,0),
                              B.kappa = B.kappa)
    nv <- mesh$n
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs)
    field  <- 1:nv
    st.volume <- diag(spde$param.inla$M0)
    expected <- c(st.volume, rep(0, nrow(locs)))
    A.pp <- rBind(Diagonal(n=nv), Ast)
    y.pp <- rep(0:1, c(nv, nrow(locs)))
    stack <- inla.stack(data=list(y=y.pp, e=expected),
                        A=list(A.pp,1),
                        effects=list(field = field, b0 = rep(1,nv+nrow(locs))))
    formula <- y ~ 0  + b0 + f(field, model=spde)
    result <- inla(as.formula(formula), family = "poisson",
                   data=inla.stack.data(stack),
                   E=inla.stack.data(stack)$e,
                   control.predictor=list(A=inla.stack.A(stack)),
                   control.inla = control.inla,
                   verbose = verbose)
    result
}
                                             
####                       
fit.ns.kappa.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["fn"]]))stop("conjectured non-stationary function must be supplied")
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    fn <- ns[["fn"]]
    B.kappa <- cbind(0,0,1,fn)
    spde <- inla.spde2.matern(mesh = mesh,
                              alpha = 2, B.tau = cbind(0,1,0,0),
                              B.kappa = B.kappa)
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    data$spde <- spde$param.inla[c("M0","M1","M2","B1","B2")]	# Encapsulation of 6 matrices
    data$area <- c(diag(data$spde$M0))
    parameters <- ns[["parameters"]]
    obj <- MakeADFun(data,parameters,random="x",DLL="nonstat",silent = verbose )
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    result <- sdreport(obj)
    result
}
                                             
####                       
fit.UN.ns.kappa.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    spde <- inla.spde2.matern(mesh = mesh,alpha = 2)
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    data$spde <- spde$param.inla[c("M0","M1","M2")]
    data$area <- c(diag(data$spde$M0))
    parameters <- ns[["parameters"]]
    obj <- MakeADFun(data,parameters,random="x",DLL="nskappaspde",silent = verbose )
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    result <- sdreport(obj)
    result
}                   

#############ns mean TMB
fit.ns.mean.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    spde <- inla.spde2.matern(mesh = mesh,alpha = 2)
    data$spde <- spde$param.inla[c("M0","M1","M2")]	
    data$area <- c(diag(data$spde$M0))
    data$dists <- as.matrix(dist(mesh$loc))
    parameters <- ns[["parameters"]] ## requires beta0 (vector), log_kappa (numeric), rho (numeric), x(vector) 
    obj <- MakeADFun(data,parameters,random="x",DLL="nonstatmean",silent = verbose )
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    result <- sdreport(obj)
    result
}


############# lgcp using TMB
fit.lgcp.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    spde <- inla.spde2.matern(mesh = mesh,alpha = 2)
    data$spde <- spde$param.inla[c("M0","M1","M2")]
    data$area <- c(diag(data$spde$M0))
    parameters <- ns[["parameters"]]
    obj <- MakeADFun(data,parameters,random="x",DLL="lgcpspde",silent = verbose )
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    result <- sdreport(obj)
    result
}



#### non-standard, non-stationary lgcp
fit.lgcp.ns <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["model"]]))stop("Please supply assumed model")
    dll.dir <- paste(system.file(package = "lgcpSPDE"), "/tmb/bin/", sep = "")
    for (i in paste(dll.dir, list.files(dll.dir), sep = "")){
        dyn.load(i)
    }
    if(ns[["model"]] == "lgcpTMB") result <- fit.lgcp.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nsmeanTMB") result <- fit.ns.mean.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nskappaTMB") result <- fit.ns.kappa.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nsUNkappaTMB") result <- fit.UN.ns.kappa.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nskappaINLA") result <- fit.ns.kappa.inla(mesh = mesh, locs = locs, ns = ns,
                                                                   control.inla = control.inla, verbose = verbose)
    result
}
