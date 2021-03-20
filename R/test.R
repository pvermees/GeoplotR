#' @title xy testing function
#' @description tests json files with bivariate decision boundaries
#' @param fname path to a json file
#' @param xlim (optional) 2-element vector with the x-axis limits
#' @param ylim (optional) 2-element vetor with the y-axis limits
#' @param log one of
#'
#' \code{''} linear axes
#' 
#' \code{'x'} semilogarithmic plot with linear y axis and logarithmic x axis
#'
#' \code{'y'} semilogarithmic plot with linear x axis and logarithmic y axis
#'
#' \code{'xy'} log-log plot
#'
#' @param short logical. If \code{TRUE}, uses short labels
#' @param polygons logical. If \code{TRUE} plots the decision
#'     polygons. If \code{TRUE}, plots the decision lines
#' @param ... additional arguments to the generic plot function
#' @return a list of text strings
#' @examples
#' fn <- system.file('TiZrY.json',package='GeoplotR')
#' xyztest(fn)
#' @export
xytest <- function(fname,xlim=NA,ylim=NA,log='',short=FALSE,polygons=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    if (any(is.na(xlim)) & any(is.na(ylim))){
        xy <- NULL
        for (pname in names(json$polygons)){
            XY <- matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            xy <-  rbind(xy,XY)
        }
        xlim <- range(xy[,1])
        ylim <- range(xy[,2])
    }
    out <- xyplot(json,xlim=xlim,ylim=ylim,show.labels=TRUE,
                  short=short,test.polygons=polygons,log=log,...)
    invisible(out)
}

#' @title xyz testing function
#' @description tests json files with ternary decision boundaries
#' @param fname path to a json file
#' @param short logical. If \code{TRUE}, uses short labels
#' @param polygons logical. If \code{TRUE} plots the decision
#'     polygons. If \code{TRUE}, plots the decision lines
#' @param ... additional arguments to the generic plot function
#' @return a list of text strings
#' @examples
#' fn <- system.file('TiZrY.json',package='GeoplotR')
#' xyztest(fn)
#' @export
xyztest <- function(fname,short=FALSE,polygons=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    out <- xyzplot(json,short=short,show.labels=TRUE,
                   test.polygons=polygons,...)
    invisible(out)
}
