#' Function that extracts the ''random fields'' of the model fitted.
#'
#' @param x A fitted model from \link{mark.pp.fit}().
#' @param mesh the mesh used in the model fit.
#' @param n.t numeric, the number of time points.
#' @param sd Logical, if \code{FALSE} means of random fields aer returned.
#' @param plot Logical, if \code{TRUE} the returned matricies (either SD or Mean of
#' random fields are plotted.
#' @param spatial.polygon Optional, if a spatial polygon of the domain is supplied, only
#' values of the random field within the domain will be returned
#' @param dims vector of length two specifying spatial resolution of projection of fields onto the mesh
#' @param ... additional graphical parameters
#'  @importFrom fields image.plot
#'
#' @export
find.fields <- function(x = NULL, mesh = NULL, n.t = NULL, sd = FALSE, plot = FALSE, spatial.polygon = NULL,
                        dims = c(100,100),...){
    if(is.null(attributes(x)$mesh) & is.null(mesh)){
        stop("no mesh has been supplied")}
    if(!is.null(attributes(x)$mesh)){mesh <- attributes(x)$mesh}else{mesh <- mesh}
    proj <- inla.mesh.projector(mesh,dims = dims)
    if(!is.null(spatial.polygon)) inside <- inwin(proj,as.owin(spatial.polygon))
    spde <-inla.spde2.matern(mesh = mesh, alpha = 2)
    fields <- summary(x)$random.names[summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy"]
    idx <- which(summary(x)$random.model=="SPDE2 model" | summary(x)$random.model=="Copy")
    n <- length(fields)
    if(!is.null(n.t)){
        t <- n.t
        means <- list()
        for (i in idx[1]:idx[n]){
            means [[i-idx[1]+1]] <- lapply(1:t, function(j) {
                r <- inla.mesh.project(proj, field = x$summary.random[[i]]$mean[1:spde$n.spde + (j-1)*spde$n.spde])
                if(!is.null(spatial.polygon)) r[!inside] <- NA; return(r)})
        }
        sds <- list()
        for (i in idx[1]:idx[n]){
            sds [[i-idx[1]+1]] <- lapply(1:t, function(j) {
                r <- inla.mesh.project(proj, field = x$summary.random[[i]]$sd[1:spde$n.spde + (j-1)*spde$n.spde])
                if(!is.null(spatial.polygon)) r[!inside] <- NA;  return(r)})
        }
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon,...)}
    }else{
        means <- list()
        for (i in idx[1]:idx[n]){
            means[[i-idx[1]+1]] <- inla.mesh.project(proj,x$summary.random[[i]]$mean)
            if(!is.null(spatial.polygon)) means[[i-idx[1]+1]][!inside] <- NA; 
        }
        sds <- list()
        for (i in idx[1]:idx[n]){
            sds[[i-idx[1]+1]] <- inla.mesh.project(proj,x$summary.random[[i]]$sd)
            if(!is.null(spatial.polygon)) sds[[i-idx[1]+1]][!inside] <- NA; 
        }
        if(plot){plot.fields( x = x, mesh = mesh, n.t = n.t, sd = sd, spatial.polygon = spatial.polygon,...)}
    }
    names(means) <- names(sds) <- fields
    ifelse(sd,return(sds),return(means))
}
