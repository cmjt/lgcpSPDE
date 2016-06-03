#' Wrapper function to create a mesh object as used in an INLA model fit
#'
#'
#' @return A ``mesh'' object used in an INLA/SPDE model
#'
#' @param locs a matrix of locations, either this or \code{spatial.polygon} must be supplied.
#' @param mesh.pars a named vertor of mesh parameters, must contain
#' \code{cutoff} length at which to cut off triangle edge lengths,
#' \code{max.edge.min} triangle edge length inside region,
#' and \code{max.edge.max} triangle edge length inside region.
#' @param spatial.polygon if supplied the spatial polygon for the domain is used to construct mesh
#' @param sphere Logical if TRUE the mesh is constructed on the unit sphere, note this is only possible if coordinates are
#' longitude and Latitude, by default FALSE
#' @param plot Logical if TRUE the triangulation is plotted
#' @export
make.mesh<-function(locs = NULL, mesh.pars = NULL, spatial.polygon = NULL, sphere = FALSE, plot = FALSE){
    # getting mesh parameters
    max.edge.min <- mesh.pars["max.edge.min"]
    max.edge.max <- mesh.pars["max.edge.max"]
    if(is.na(max.edge.max)){max.edge <- max.edge.min}else{max.edge <- c(max.edge.min,max.edge.max)}
    cutoff <- mesh.pars["cutoff"]
    mesh.pars<-c(max.edge,cutoff)
    if(!sphere){
        if(!is.null(spatial.polygon)){
            # creates triangulation based on a spatial polygon of the domain
            boundary <- inla.sp2segment(spatial.polygon)
            mesh <- inla.mesh.2d(boundary = boundary, max.edge = max.edge, cutoff = cutoff)
        } else {
            loc <- locs
            # creates triangulation based on the locations of the point pattern
            mesh <- inla.mesh.2d(loc = locs, max.edge = max.edge, cutoff = cutoff)
        }}
    if(sphere){
        if(!is.null(spatial.polygon)){
            # creates triangulation based on a spatial polygon of the domain projected onto a sphere
            boundary <- inla.sp2segment(spatial.polygon)
            boundary$loc <- inla.mesh.map(boundary$loc, projection="longlat", inverse=TRUE)
            mesh <- inla.mesh.2d(boundary = boundary, max.edge = max.edge, cutoff = cutoff)
        } else {
            # creates triangulation based on the locations of the point pattern projected onto a sphere
            locs <- inla.mesh.map(locs, projection="longlat", inverse=TRUE)
            mesh <- inla.mesh.2d(loc = locs, max.edge = max.edge, cutoff = cutoff)
        }}
    if(plot) plot.mesh(mesh)
    mesh
}
