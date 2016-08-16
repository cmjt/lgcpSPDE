#' Function to simulate from a spatial/spatio temporal LGCP (SPDE model)
#'
#'
#' @return A named matrix (or a list of matricies if spatio-temporal) of point locations and
#' (if a marked point pattern is simulated) a mark values
#' @param spatial.polygon the spatial polygon for the domain is used to construct the delauney traingulation
#' @param mesh.pars a named vertor of mesh parameters, must contain
#' \code{cutoff} length at which to cut off triangle edge lengths,
#' \code{min} triangle edge length inside region,
#' and \code{max} triangle edge length inside region.
#' @param mu numeric, the intercept term to simulate a LGCP, by default is 0.
#' @param kappa a numeric constant, parameter of the SPDE model.
#' @param sigma2 a numeric constant, parameter of the SPDE model, by default this is 0.05.
#' @param n a numeric constant defining the number of time points, by default 1.
#' @param rho the ar1 correlation coefficient for spatio-temporal samples,
#' by default this is 0.9.
#' @param mark Logical, if TRUE a marked point pattern is simulated
#' @param beta a scalar, this, the interaction parameter describing the dependance between the mark and point locations 
#' @param mark.function a function of 2D spatial coordinates which describes the spatial process
#' specific to the mark, by default this is \code{function(x,y) cos(x) - sin(y)}.
#' @param seed seed for the simulation, by default this is 1
#' @export

rlgcpspde<-function (spatial.polygon = NULL, mesh.pars = NULL, mu = 0, kappa = NULL, sigma2 = 0.05 , n = 1, rho = 0.9, mark = FALSE, beta = NULL, mark.function = function(x,y) cos(x) - sin(y), seed = 1){
    mesh <- make.mesh(mesh.pars = mesh.pars, spatial.polygon = spatial.polygon)
    locs <- mesh$loc
    sample <- mu + rgeospde(locs = locs, mesh = mesh, kappa = kappa, sigma2 = sigma2, n = n, rho = rho, seed = seed)
    proj<-inla.mesh.projector(mesh = mesh)
    w <- as.owin(spatial.polygon)
    y0 <- x0 <- seq(w$xrange[1], w$xrange[2],length=length(proj$x))
    if(mark){
        mark.im <- outer(x0,y0, mark.function)
        if(ncol(sample)==1){
            logLambda <- matrix(inla.mesh.project(proj, sample),length(proj$x),length(proj$y))
            set.seed(seed)
            pp <-rpoispp(as.im(exp(logLambda), W = w))[w]
            loc <- cbind(pp$x,pp$y)
            mark <- mark.im[Reduce('cbind', nearest.pixel(loc[,1],loc[,2],im(mark.im, x0, y0)))]
            ppstruct <- logLambda[Reduce('cbind', nearest.pixel(loc[,1],loc[,2],im(logLambda, x0, y0)))]
            mark <- mark + beta*ppstruct
            result <- cbind(x = loc[,1],y = loc[,2], mark = mark)
        }else{
            logLambda <- lapply(1:ncol(sample), function(j) {
                r <- matrix(inla.mesh.project(proj, sample[,j]),length(proj$x),length(proj$y))
                return(r)})
            pp <- lapply(1:ncol(sample), function(j) {
                set.seed(seed)
                r <-rpoispp(as.im(exp(logLambda[[j]]), W = w))[w]
                return(r)})
            loc <- lapply(pp, function(x) cbind(x$x,x$y))
            ppstruct <- sapply(1:length(pp), function(i) {
                logLambda[[i]][Reduce('cbind', nearest.pixel(loc[[i]][,1],loc[[i]][,2],im(logLambda[[i]], x0, y0)))]})
            mark <- sapply(1:length(pp), function(i) {
                mark.im[Reduce('cbind', nearest.pixel(loc[[i]][,1],loc[[i]][,2],im(mark.im, x0, y0)))]})
            mark <-  sapply(1:length(pp),function(i) mark[[i]] + beta*ppstruct[[i]])
            result <- sapply(1:length(pp), function(i) cbind(x = loc[[i]][,1],y = loc[[i]][,2], mark = mark[[i]]))
            }}else{
                 if(ncol(sample)==1){
                     Lambda <- as.im(exp(matrix(inla.mesh.project(proj, sample),length(proj$x),length(proj$y))),W = w)
                     set.seed(seed)
                     pp <-rpoispp(Lambda)[w]
                     loc <- cbind(pp$x,pp$y)
                     result <- cbind(x = loc[,1],y = loc[,2])
                 }else{
                     logLambda <- lapply(1:ncol(sample), function(j) {
                         r <- matrix(inla.mesh.project(proj, sample[,j]),length(proj$x),length(proj$y))
                         return(r)})
                     pp <- lapply(1:ncol(sample), function(j) {
                         set.seed(seed)
                         r <-rpoispp(as.im(exp(logLambda[[j]]), W = w))[w]
                         return(r)})
                     loc <- lapply(pp, function(x) cbind(x$x,x$y))
                     result <- sapply(1:length(pp), function(i) cbind(x = loc[[i]][,1],y = loc[[i]][,2]))
                 }
            }
    return(result)
}


            
            



