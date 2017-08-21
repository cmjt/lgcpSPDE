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
    if(length(resp.c)==0)stop("mesh needs to be constructed using point pattern locations")
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
fit.oscillate.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    spde <- inla.spde2.matern(mesh = mesh,
                              alpha = 2)
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    if(length(resp.c)==0)stop("mesh needs to be constructed using point pattern locations")
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    data$spde <- spde$param.inla[c("M0","M1","M2")]	# Encapsulation of 6 matrices
    data$area <- c(diag(data$spde$M0))
    parameters <- ns[["parameters"]]
    obj <- MakeADFun(data,parameters,random="x",DLL="osns",silent = verbose )
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    result <- sdreport(obj)
    result
}
                                             
####                       
fit.UN.ns.kappa.TMB <- function(mesh = NULL, locs = NULL, ns = NULL, control.inla = NULL, verbose = NULL){
    if(is.null(ns[["parameters"]]))stop("TMB requires parameter starting values")
    spde <- inla.spde2.matern(mesh = mesh,alpha = 2,B.tau = cbind(0,1,0,0),
                              B.kappa = cbind(0,0,1,1))
    resp <- numeric(mesh$n)
    resp.c <-as.vector(table(mesh$idx$loc))
    if(length(resp.c)==0)stop("mesh needs to be constructed using point pattern locations")
    resp[unique(mesh$idx$loc)]<- resp.c
    meshidxloc <- 1:mesh$n
    data <- list(resp = resp, meshidxloc=meshidxloc)
    data$spde <- spde$param.inla[c("M0","M1","M2","B2")]
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
    if(length(resp.c)==0)stop("mesh needs to be constructed using point pattern locations")
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
    if(length(resp.c)==0)stop("mesh needs to be constructed using point pattern locations")
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
    ##if(ns[["model"]] == "nsmeanTMB") result <- fit.ns.mean.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nskappaTMB") result <- fit.ns.kappa.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nsUNkappaTMB") result <- fit.UN.ns.kappa.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "oscillateTMB") result <- fit.oscillate.TMB(mesh = mesh, locs = locs, ns = ns, verbose = !verbose)
    if(ns[["model"]] == "nskappaINLA") result <- fit.ns.kappa.inla(mesh = mesh, locs = locs, ns = ns,
                                                                   control.inla = control.inla, verbose = verbose)
    result
}

#### Haakon's code to use pc priors for parameters of latent field

inla.spde2.matern.new = function(mesh, alpha=2, prior.pc.rho, prior.pc.sig){
    d = INLA:::inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
    nu = alpha-d/2
    kappa0 = log(8*nu)/2
    tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
    spde   = inla.spde2.matern(mesh = mesh,
                               B.tau   = cbind(tau0,   nu,  -1),
                               B.kappa = cbind(kappa0, -1, 0))
    param = c(prior.pc.rho, prior.pc.sig)
    spde$f$hyper.default$theta1$prior = "pcspdega"
    spde$f$hyper.default$theta1$param = param
    spde$f$hyper.default$theta1$initial = log(prior.pc.rho[1])+1
    spde$f$hyper.default$theta2$initial = log(prior.pc.sig[1])-1
    
  # End and return
    return(invisible(spde))  
}


### brief sim.fun for bird paper in based on Chapter 3 of spde tutorial
## (http://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/spde-tutorial.pdf)

sim.fun.bird.paper <- function(n.sim = 1,n = NULL,alpha = NULL,m.var = NULL, kappa = NULL,beta = NULL,e.sd = NULL,
                               b.param = c(0,10)){
    alpha = alpha
    m.var = m.var
    kappa = kappa
    beta = beta
    n = n
    e.sd = e.sd
    locs = cbind(runif(n), runif(n))
    fix1 =fix2 = bet = range1 = sd1 = range2 = sd2 = matrix(NA,nrow = n.sim,ncol = 6)
    colnames(fix1) = colnames(fix2) = colnames(bet)= c("mean","sd","0.25 quant","0.5 quant","0.975quant","mode")
    colnames(range1)=colnames(range2)=colnames(sd1)=colnames(sd2) =  c("mean","sd","0.25 quant","0.5 quant","0.975quant","mode")
    for(i in 1:n.sim){
        z1 = rMatern(1, locs, kappa[1], m.var[1])
        z2 = rMatern(1, locs, kappa[2], m.var[2])
                                        #obvs samples
        
        y1 = exp(alpha[1] + z1[1:n] + rnorm(n, 0, e.sd[1]))
        y2 = alpha[2] + beta[1] * z1[1:n] + z2[1:n] + rnorm(n, 0, e.sd[2])
        y2 = ifelse(y2>=0,1,0)
        mesh <- inla.mesh.2d(locs, max.edge=c(0.05, 0.2),
                             offset=c(0.05, 0.3), cutoff=0.01)
        spde <- inla.spde2.pcmatern(
            mesh=mesh, alpha=2, ### mesh and smoothness parameter
            prior.range=c(0.05, 0.01), ### P(practic.range<0.05)=0.01
            prior.sigma=c(1, 0.01)) ### P(sigma>1)=0.01
        hc1 <- list(theta=list(prior='normal', param = b.param))
        form <- y ~ 0 + intercept1 + intercept2 +
            f(s1, model=spde) + f(s2, model=spde) +  f(s12, copy="s1", fixed=FALSE, hyper=hc1)
        A <- inla.spde.make.A(mesh, locs)
        stack1 <- inla.stack(
            data=list(y=cbind(as.vector(y1), NA)), A=list(A),
            effects=list(list(intercept1=1, s1=1:spde$n.spde)))
        stack2 <- inla.stack(
            data=list(y=cbind(NA, as.vector(y2))), A=list(A),
            effects=list(list(intercept2=1, s2=1:spde$n.spde,
                              s12=1:spde$n.spde)))
        stack <- inla.stack(stack1, stack2)
        #theta.ini = c(log(1/e.sd[1]^2),c(log(sqrt(8)/kappa), log(sqrt(m.var)))[c(1,3,2,4)], beta)
        result <- inla(form, c("gamma","binomial"),
                       data=inla.stack.data(stack),
                       control.predictor=list(A=inla.stack.A(stack)),
                       #control.mode=list(theta=theta.ini, restart=TRUE),
                       control.inla=list(int.strategy='eb'))
        fix1[i,] = summary(result)$fixed[1,1:6]
        fix2[i,] = summary(result)$fixed[2,1:6]
        bet[i,] = as.matrix(summary(result)$hyperpar)[6,1:6]
        range1[i,] = as.matrix(summary(result)$hyperpar)[2,1:6]
        range2[i,] = as.matrix(summary(result)$hyperpar)[4,1:6]
        sd1[i,] = as.matrix(summary(result)$hyperpar)[3,1:6]
        sd2[i,] = as.matrix(summary(result)$hyperpar)[5,1:6]
        cat("iteration",i,"\n")
        }
    return(list(intercpt1 = fix1, intercept2 = fix2, Beta = bet,
                range1 = range1, range2 = range2,sd1 = sd1, sd2 = sd2))
     ##return(result)
}

## function from http://www.math.ntnu.no/inla/r-inla.org/tutorials/spde/R/spde-tutorial-functions.R

rMatern <- function(n, coords, kappa, variance, nu=1) {
    m <- as.matrix(dist(coords))
    m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
    diag(m) <- 1
    return(drop(crossprod(chol(variance*m),
                          matrix(rnorm(nrow(coords)*n), ncol=n))))
}
