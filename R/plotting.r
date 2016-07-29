##plots fields eith spatio-temporal or spatial; function that calls a function that calls it!!
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

## plot mesh
plot.mesh <- function(x){
    plot(x,main="",asp=1,draw.segment = FALSE)
    if (!is.null(x$segm$bnd)) 
                lines(x$segm$bnd, x$loc, lwd = 2,col = 1)
    if (!is.null(x$segm$int)) 
                lines(x$segm$int, x$loc, lwd = 2,col = 1)
}

