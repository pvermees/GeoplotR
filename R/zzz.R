.onLoad <- function(libname, pkgname){
    utils::data(training,test,package=pkgname,envir=parent.env(environment()))
}
