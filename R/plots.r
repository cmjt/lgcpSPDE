#' Plotting the estimated means of the random fields or parameter densities
#'
#' Plots the estimated random fields or parameter posterior densities from an object returned by
#' \link{mark.pp.fit}().
#'
#' @param x A fitted model from \link{mark.pp.fit}().
#' @param mesh the mesh used in the model fit.
#' densities for the model parameters are plotted.
#' @importFrom fields image.plot
#' 
#'
#' @export
plot.fields <- function(x = NULL, mesh = NULL, n.t = NULL, sd = FALSE){
    proj <- inla.mesh.projector(mesh)
    fields <- names(x$summary.random)
    n <- length(fields)
    par(mfrow=c(3,3))
    if(!is.null(n.t)){
        rfs <- find.fields(x = x, mesh = mesh, n.t = n.t, sd = sd)
        t <- n.t
        for(i in 1:n){
            for(j in 1:t){ image.plot(proj$x,proj$y,rfs[[i]][[j]],axes=FALSE,xlab="",ylab="", main = paste(fields[i], "time", j,  sep = " "))}
        }
    }else{
        rfs <- find.fields(x = x, mesh = mesh, sd = sd)
        for(i in 1:n){
            image.plot(proj$x,proj$y,rfs[[i]],axes=FALSE,xlab="",ylab="", main = fields[i])
        }
    }
}
