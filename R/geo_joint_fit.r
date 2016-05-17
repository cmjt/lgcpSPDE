#' Function to fit a joint spatio-temporal model to geo-statistical data where one spatio-
#' temporal component is shared between the responses
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs a matrix of observation locations, where each row corresponds to the observation. 
#' @param response a matrix of response variables, each row corresponds to the spatial locations
#' in \code{locs}, the first column is the first response as a result of a
#' spatio-temporal process we want to ``copy'' into the second column response.
#' @param temp a numeric vector specifying a temporal index for each observation (starting at 1.....T)
#' @param family a character vector specifying the assumed likelihoods of the response, by default is c("gaussian,"gaussian").
#' @param prior.rho prior for the temporal correlation coefficient, by default a \code{pcprior} is used with \code{param=c(0-0.9)}. 
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
geo.joint.fit <- function(mesh = NULL,  locs = NULL, response = NULL, temp = NULL, family = c("gaussian","gaussian"), prior.rho = list(theta = list(prior='pccor1', param = c(0, 0.9))), hyper = list(theta=list(prior='normal', param=c(0,10))), verbose = FALSE){
    response.1 <- response[,1]
    response.2 <- response[,2]
    # spde model for the spatial random field
    spde <- inla.spde2.matern(mesh, alpha = 2)
    if(!is.null(temp)){
        temp <- temp # temporal dimention
        k <- max(temp)
        field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde, n.group = k)
        field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde, n.group = k)
        copy.field.1 <- inla.spde.make.index('copy.field.1', n.spde = spde$n.spde, n.group = k)
        # make the projector matrix (as both responses are observed at the same locations
        # this is the same for both
        A <- inla.spde.make.A(mesh = mesh, loc = locs, n.group = k, group = temp)
        #create data stacks
        stack.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                              A = list(A), 
                              tag = 'response.1',
                              effects = list(field.1 = field.1))
        stack.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                              A = list(A,A), 
                              tag = 'response.2',
                              effects = list(field.2 = field.2, copy.field.1 = copy.field.1))
        stack<-inla.stack(stack.1,stack.2)
        # temporal model and priors
        control <- list(model = 'ar1', hyper = prior.rho)
        #formula and call to inla
        formula <- y ~ 0 + f(field.1,model = spde, group = field.1.group, control.group = control) +
            f(field.2, model = spde, group = field.2.group, control.group = control) +
            f(copy.field.1, copy = 'field.1',fixed = FALSE, hyper = hyper)
    }else{
        field.1 <- inla.spde.make.index('field.1', n.spde = spde$n.spde)
        field.2 <- inla.spde.make.index('field.2', n.spde = spde$n.spde)
        copy.field.1 <- inla.spde.make.index('copy.field.1', n.spde = spde$n.spde)
        # make the projector matrix (as both responses are observed at the same locations
        # this is the same for both
        A <- inla.spde.make.A(mesh = mesh, loc = locs)
        #create data stacks
        stack.1 <- inla.stack(data=list(y=cbind(response.1,NA)),
                              A = list(A), 
                              tag = 'response.1',
                              effects = list(field.1 = field.1))
        stack.2 <- inla.stack(data=list(y=cbind(NA,response.2)),
                              A = list(A,A), 
                              tag = 'response.2',
                              effects = list(field.2 = field.2, copy.field.1 = copy.field.1))
        stack<-inla.stack(stack.1,stack.2)
        #formula and call to inla
        formula <- y ~ 0 + f(field.1,model = spde) +
            f(field.2, model = spde) +
            f(copy.field.1, copy = 'field.1',fixed = FALSE, hyper = hyper)
    }
        
    result = inla(formula, data=inla.stack.data(stack),
                  family = family,only.hyperparam = FALSE,
                  control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
                  control.compute = list(config = TRUE),
                  verbose = verbose)
    attributes <- list()
    attributes$stack <- stack
    attributes(result) <- c(attributes(result), attributes)
    result
}
