#' @title xy testing function
#' @description tests json files with bivariate decision boundaries
#' @param fname path to a json file
#' @param xlim (optional) 2-element vector with the x-axis limits
#' @param ylim (optional) 2-element vetor with the y-axis limits
#' @param xlab (optional) x-axis label
#' @param ylab (optional) y-axis label
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
#' @param show.labels logical. If \code{TRUE}, labels the decision
#'     fields.
#' @param short logical. If \code{TRUE}, uses short labels
#' @param polygons logical. If \code{TRUE} plots the decision
#'     polygons. If \code{TRUE}, plots the decision lines
#' @param smooth logical. If \code{TRUE}, plots lines as
#'     \code{xspline}s.
#' @param ... additional arguments to the generic plot function
#' @return a list of text strings
#' @examples
#' fn <- system.file('TiZrY.json',package='GeoplotR')
#' xyztest(fn)
#' @export
xytest <- function(fname,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,log='',
                   show.labels=TRUE,short=FALSE,smooth=FALSE,polygons=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    if (is.null(xlim) | is.null(ylim)){
        xy <- NULL
        for (pname in names(json$lines)){
            XY <- matrix(unlist(json$lines[[pname]]),ncol=2,byrow=TRUE)
            xy <-  rbind(xy,XY)
        }
        if (is.null(xlim)) xlim <- range(xy[,1])
        if (is.null(ylim)) ylim <- range(xy[,2])
    }
    out <- xyplot(json,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,
                  show.labels=show.labels,short=short,
                  test.polygons=polygons,smooth=smooth,log=log,...)
    graphics::title(main=fname)
    invisible(out)
}

#' @title xyz testing function
#' @description tests json files with ternary decision boundaries
#' @param fname path to a json file
#' @param xyzlab (optional) 3-element vector of corner labels.
#' @param show.labels logical. If \code{TRUE}, labels the decision
#'     fields.
#' @param short logical. If \code{TRUE}, uses short labels.
#' @param polygons logical. If \code{TRUE} plots the decision
#'     polygons. If \code{TRUE}, plots the decision lines
#' @param smooth logical. If \code{TRUE}, plots lines as
#'     \code{xspline}s.
#' @param ... additional arguments to the generic plot function
#' @return a list of text strings
#' @examples
#' fn <- system.file('TiZrY.json',package='GeoplotR')
#' xyztest(fn)
#' @export
xyztest <- function(fname,xyzlab=NULL,show.labels=TRUE,short=FALSE,
                    polygons=FALSE,smooth=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    out <- xyzplot(json,xyzlab=xyzlab,show.labels=show.labels,
                   short=short,test.polygons=polygons,smooth=smooth,...)
    graphics::title(main=fname)
    invisible(out)
}
