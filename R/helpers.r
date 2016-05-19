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

##find fields eith spatio temporal or spatial
find.fields <- function(x = NULL, mesh = NULL, n.t = NULL, sd = FALSE){
    proj <- inla.mesh.projector(mesh)
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    fields <- names(x$summary.random)
    n <- length(fields)
    if(!is.null(n.t)){
        t <- n.t
        means <- list()
        for (i in 1:n){
            means [[i]] <- lapply(1:t, function(j) { r <- inla.mesh.project(proj, field = x$summary.random[[i]]$mean[1:spde$n.spde + (j-1)*spde$n.spde]); return(r)})
        }
        sds <- list()
        for (i in 1:n){
            sds [[i]] <- lapply(1:t, function(j) {r <- inla.mesh.project(proj, field = x$summary.random[[i]]$sd[1:spde$n.spde + (j-1)*spde$n.spde]); return(r)})
        }
    }else{
        means <- list()
        for (i in 1:n){
            means[[i]] <- inla.mesh.project(proj,x$summary.random[[i]]$mean)
            }
        sds <- list()
        for (i in 1:n){
            sds[[i]] <- inla.mesh.project(proj,x$summary.random[[i]]$sd)
            }
    }
    ifelse(sd,return(sds),return(means))
}

## model fit for ant model (will generalise at some point) i.e. spatio temporal with two marks that is, p/a or ant nests, and p/a of pests

ants.pp.fit <- function(mesh = NULL, locs=NULL, t.index = NULL, marks = NULL, covariates = NULL, mark.family = c("binomial","binomial"), verbose = FALSE,
                        prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))),hyper = list(theta=list(prior='normal', param=c(0,10)))){
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    # number of observations
    n <- nrow(locs)
    # number of mesh nodes
    nv <- mesh$n
    temp <- t.index # temporal dimention
    k <- (mesh.t <- inla.mesh.1d(temp))$n # number of groups
    ## the response for the point pattern locations
    y.pp <- rep(0:1, c(k * nv, n))
    mark.1 <- marks[,1]
    mark.2 <- marks[,2]
    ## create projection matrix for loacations
    Ast <- inla.spde.make.A(mesh = mesh, loc = locs, group = temp, n.group = k)
    ##effect for LGCP used for point pattern
    st.volume <- diag(kronecker(Diagonal(n = k),spde$param.inla$M0))
    expected <- c(st.volume, rep(0, n))
    field.pp <- inla.spde.make.index('field.pp', n.spde = spde$n.spde, group = temp, n.group = k)
    field.mark.1 <- inla.spde.make.index('field.mark.1', n.spde = spde$n.spde, group = temp, n.group = k)
    copy.field.mark.1 <- inla.spde.make.index('copy.field.mark.1', n.spde = spde$n.spde, group = temp, n.group = k)
    copy.field.1 <- inla.spde.make.index('copy.field.1', n.spde = spde$n.spde, group = temp, n.group = k)
    copy.field.2 <- inla.spde.make.index('copy.field.2', n.spde = spde$n.spde, group = temp, n.group = k)
                                        # temporal model "ar1"
    ctr.g <- list(model='ar1',param = prior.rho)
    stk.pp <- inla.stack(data=list(y=cbind(y.pp,NA,NA), e=expected),
                         A=list(rBind(Diagonal(n=k*nv), Ast)),
                         effects=list(field.pp = field.pp))
    stk.mark.1 <- inla.stack(data=list(y=cbind(NA,mark.1,NA)),
                             A=list(Ast, Ast),
                             effects=list(field.mark.1 = field.mark.1, copy.field.1 = copy.field.1))
    stk.mark.2 <- inla.stack(data=list(y=cbind(NA,NA,mark.2)),
                             A=list(Ast, Ast ),
                             effects=list(copy.field.mark.1 = copy.field.mark.1, copy.field.2 = copy.field.2))
    stack <- inla.stack(stk.pp,stk.mark.1, stk.mark.2)
    formula <- y ~ 0 + f(field.pp, model=spde, group = field.pp.group, control.group=ctr.g) +
        f(field.mark.1, model=spde, group = field.mark.1.group , control.group=ctr.g) +
        f(copy.field.1, copy = "field.pp", fixed=FALSE, hyper = hyper ) +
        f(copy.field.2, copy = "field.pp", fixed=FALSE, hyper = hyper ) +
        f(copy.field.mark.1, copy = "field.mark.1", fixed=FALSE, hyper = hyper )
    ##call to inla
    result <- inla(formula, family = c("poisson",mark.family),
            data=inla.stack.data(stack),
            E=inla.stack.data(stack)$e,
            control.predictor=list(A=inla.stack.A(stack)),
            control.inla=list(strategy='gaussian',int.strategy = 'eb'),
            verbose = verbose)
    
    result
}
                                      
    

                                            
    
