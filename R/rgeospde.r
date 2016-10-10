#' Function to simulate from a spatial/spatio temporal SPDE model
#'
#'
#' @return A matrix of values each column a set of observations at each time point
#'
#' @param locs a matrix of locations at which the values are to be simulated
#' @param mesh a ``mesh'' object i.e. delauney triangulation of the domain, an
#' object returned by \link{make.mesh}.
#' @param kappa a numeric constant, parameter of the SPDE model.
#' @param sigma2 a numeric constant, parameter of the SPDE model, by default this is 1.
#' @param n a numeric constant defining the number of time points, by default 1.
#' @param rho the ar1 correlation coefficient for spatio-temporal samples,
#' by default this is 0.9.
#' @param seed seed for the simulation, by default this is 1
#' @param non.stat a named list the first element \code{fn} the spatial function which kappa varies with
#' and \code{theta} a vctor of length 3 specifying the theta values of the non-stationary model
#' @export
rgeospde <- function( locs = NULL,mesh = NULL, kappa = NULL, sigma2 = 1,n = 1, rho = 0.9, seed = 1, non.stat = NULL){
    if(!is.null(non.stat) & is.null(non.stat[["oscillate"]])){
        fn <- non.stat[["fn"]]
        theta <- non.stat[["theta"]]
        B.kappa <- cbind(0,0,1,fn)
        spde <- inla.spde2.matern(mesh = mesh,
                                 alpha = 2, B.tau = cbind(0,1,0,0),
                                 B.kappa = B.kappa)
    }else{
        spde <- inla.spde2.matern(mesh = mesh, alpha=2)
        theta <- c(-0.5 * log(4 * pi * sigma2 * kappa^2), log(kappa))}
    if(!is.null(non.stat[["oscillate"]])){
        thetaOSC <- non.stat[["oscillate"]]
        Q <- kappa^4 *spde$param.inla$M0 +
                       2*kappa^2*cos(pi*thetaOSC)*spde$param.inla$M1 +
                                   spde$param.inla$M2
    }else{
        Q <- inla.spde2.precision(spde, theta)
    } # PRECISION MATRIX
    #projector matrix
    A <- inla.mesh.project(mesh = mesh,loc = locs)$A
    if (n == 1) {
        # sample at the mesh nodes
        x <- inla.qsample(Q = Q, seed = seed, constr = spde$f$extraconstr)
        x <- drop(A%*%x)
    } else {
        x.t <- inla.qsample (n = n, Q = Q, seed = seed, constr = spde$f$extraconstr)
        x.t <- drop(A%*%x.t)
        x <- x.t
        for (j in 2:n){
            x[,j] <- rho*x[,j-1] + sqrt(1-rho^2)*x.t[,j]
        }
    }
    x
}
