#' Function that extracts the ''random fields'' of the model fitted.
#'
#' Plots the estimated random fields or parameter posterior densities from an object returned by
#' \link{mark.pp.fit}().
#'
#' @param x A fitted model from \link{mark.pp.fit}().
#' @param mesh the mesh used in the model fit.
#' @param n.t numeric, the number of time points.
#' @param sd Logical, if \code{FALSE} means of random fields aer returned.
#' @param plot Logical, if \code{TRUE} the returned matricies (either SD or Mean of
#' random fields are plotted.
#' @importFrom spatstat as.owin
#'  @importFrom fields image.plot
#'
#' @export
find.fields <- function(x = NULL, mesh = NULL, n.t = NULL, sd = FALSE, plot = FALSE, spatial.polygon = NULL){
    if(is.null(attributes(x)$mesh) & is.null(mesh)){
        stop("no mesh has been supplied")}
    if(!is.null(attributes(x)$mesh)){mesh <- attributes(x)$mesh}else{mesh <- mesh}
    proj <- inla.mesh.projector(mesh)
    if(!is.null(spatial.polygon)) inside <- inwin(proj,as.owin(spatial.polygon))
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    fields <- names(x$summary.random)
    n <- length(fields)
    if(!is.null(n.t)){
        t <- n.t
        means <- list()
        for (i in 1:n){
            means [[i]] <- lapply(1:t, function(j) { r <- inla.mesh.project(proj, field = x$summary.random[[i]]$mean[1:spde$n.spde + (j-1)*spde$n.spde]);  if(!is.null(spatial.polygon)) r[!inside] <- NA; return(r)})
        }
        sds <- list()
        for (i in 1:n){
            sds [[i]] <- lapply(1:t, function(j) {r <- inla.mesh.project(proj, field = x$summary.random[[i]]$sd[1:spde$n.spde + (j-1)*spde$n.spde]);if(!is.null(spatial.polygon)) r[!inside] <- NA;  return(r)})
        }
        if(!is.null(spatial.polygon)) for(i in 1:n){sds[[i]][!inside] <- NA}
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon)}
    }else{
        means <- list()
        for (i in 1:n){
            means[[i]] <- inla.mesh.project(proj,x$summary.random[[i]]$mean)
            if(!is.null(spatial.polygon)) means[[i]][!inside] <- NA; 
            }
        sds <- list()
        for (i in 1:n){
            sds[[i]] <- inla.mesh.project(proj,x$summary.random[[i]]$sd)
            if(!is.null(spatial.polygon)) sds[[i]][!inside] <- NA; 
            }
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon)}
    }
    ifelse(sd,return(sds),return(means))
}
