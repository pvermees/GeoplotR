#' @title Hollocher et al. (2012)
#' @name Hollocher2012
#' @rdname Hollocher2012
#' @description Discriminants for MORB, ocean island, and a variety of
#'     arc-type basalts
#' @param Nb vector with Nb concentrations (ppm)
#' @param La vector with La concentrations (ppm)
#' @param Yb vector with Yb concentrations (ppm)
#' @param Th vector with Th concentrations (ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     fields with their respective tectonic affinities.
#' @param short logical. If \code{TRUE}, uses abbreviated labels.
#' @param ... additional arguments to the generic \code{points}
#'     function.
#' @references Hollocher et al., 2012 Discriminants for MORB, ocean
#'     island, and a variety of arc-type basalts
#' @examples
#' data(test,package='GeoplotR')
#' @return a vector with tectonic affinities
NULL

#' @rdname Hollocher2012
#' @examples
#' NbLaYb(Nb=test[,'Nb'],La=test[,'La'],Yb=test[,'Yb'])
#' @export
NbLaYb <- function(Nb=NULL,La=NULL,Yb=NULL,xlim=NULL,ylim=NULL,
                   show.labels=TRUE,short=FALSE,...){
    if (is.null(xlim)) xlim <- getlimits(x=La/Yb,m=0.5,M=60)
    if (is.null(ylim)) ylim <- getlimits(x=Nb/La,m=0.1,M=2.5)
    invisible(xyplot(json=.NbLaYb,X=La/Yb,Y=Nb/La,log='xy',
                     xlim=xlim,ylim=ylim,
                     show.labels=show.labels,short=short,...))
}

#' @rdname Hollocher2012
#' @examples
#' ThNbLaYb(Th=test[,'Th'],Nb=test[,'Nb'],La=test[,'La'],Yb=test[,'Yb'])
#' @export
ThNbLaYb <- function(Th=NULL,Nb=NULL,La=NULL,Yb=NULL,
                   xlim=NULL,ylim=NULL,show.labels=TRUE,short=FALSE,...){
    if (is.null(xlim)) xlim <- getlimits(x=La/Yb,m=0.4,M=80)
    if (is.null(ylim)) ylim <- getlimits(x=Th/Nb,m=0.02,M=5)
    invisible(xyplot(json=.ThNbLaYb,X=La/Yb,Y=Th/Nb,log='xy',
                     xlim=xlim,ylim=ylim,show.labels=show.labels,short=short,...))
}
