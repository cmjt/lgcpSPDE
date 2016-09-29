#' Compiling TMB C++ templates
#'
#' Compiles the \code{lgcpSPDE} TMB templates into a shared object
#' file. This must be done a single time following installation or
#' updating of the package.
#'
#' @export
compile.lgcpSPDE <- function(){
    wd <- getwd()
    dir <- paste(system.file(package = "lgcpSPDE"), "/tmb/src", sep = "")
    setwd(dir)
    if (!dir.exists("../bin")){
        dir.create("../bin")
    }
    files <- strsplit(list.files(), "[.]")
    base <- sapply(files, function(x) x[1])
    ext <- sapply(files, function(x) x[2]) 
    for (i in base[ext == "cpp"]){
        compile(paste(i, ".cpp", sep = ""))
        unlink(paste(i, ".o", sep = ""))
        file.rename(paste(i, ".so", sep = ""),
                    paste("../bin/", i, ".so", sep = ""))
    }
    setwd(wd)
}

#' @import TMB
