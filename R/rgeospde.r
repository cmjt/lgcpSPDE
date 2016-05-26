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
#' @export
rgeospde <- function( locs = NULL,mesh = NULL, kappa = NULL, sigma2 = 1,n = 1, rho = 0.9, seed = 1){
    # define the spde model for the domain
    spde <- inla.spde2.matern(mesh = mesh, alpha=2)
    # define theta
    theta <- c(-0.5 * log(4 * pi * sigma2 * kappa^2), log(kappa))
    # precision matrix
    Q <- inla.spde2.precision(spde, theta)
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
