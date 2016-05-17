#' Function to fit a intercept only spatial geo-statistical ''prediction' model
#'
#'
#' @return A \code{inla} result object
#'
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param locs.o a matrix of observation locations, where each row corresponds to the observation.
#' @param locs.p a matrix of observation locations, where each row corresponds to the prediction location. 
#' @param response a vector of response variables, each corresponds to the spatial locations
#' in \code{locs}.
#' @param family a character specifying the assumed likelihood of the response, by default is Gaussian.
#' @param verbose Logical if \code{TRUE} model fit is output to screen.
geo.fit <- function(mesh = NULL,  locs.o = NULL, locs.p = NULL, response = NULL, family = "gaussian",verbose = FALSE){
    # spde model for the spatial random field
    spde <- inla.spde2.matern(mesh, alpha = 2)
    # index and data stack for observations
    index.o <- inla.spde.make.A(mesh = mesh, loc = locs.o)
    index.p <- inla.spde.make.A(mesh = mesh, loc = locs.p)
    stk.obvs <- inla.stack(data=list(y = response),
                        A=list(index.o), tag='observation',
                        effects=list(field = 1:mesh$n))
    # then for the prediction
    stk.prd <- inla.stack(data=list(y=NA), A=list(index.p),
                          effects=list(field = 1:spde$n.spde), tag='prediction')
    # combine data stacks
    stack<-inla.stack(stk.obvs,stk.prd)
    result <- inla(y ~ 0 + f(field,model = spde), family = family ,data = inla.stack.data(stack),
                   control.predictor=list(compute=TRUE, A=inla.stack.A(stack)),
                   verbose = verbose)
    attributes <- list()
    attributes$stack <- stack
    attributes(result) <- c(attributes(result), attributes)
    result
}
