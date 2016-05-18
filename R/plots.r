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
plot.fields <- function(x, mesh,  plot.parameters){
    proj <- inla.mesh.projector(mesh)
    fields <- names(x$summary.random)
    n <- length(fields)
    par(mfrow=c(n,n))
    for(i in 1:n){
        image.plot(proj$x,proj$y,inla.mesh.project(proj,x$summary.random[[i]]$mean),axes=FALSE,xlab="",ylab="", main = paste(fields[i], "mean", sep = " "))
        }
    for(i in 1:n){
        image.plot(proj$x,proj$y,inla.mesh.project(proj,x$summary.random[[i]]$sd),axes=FALSE,xlab="",ylab="", main = paste(fields[i], "SD", sep = " "))
        }
    
}
