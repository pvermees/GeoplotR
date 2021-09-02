#' @title Pearce et al. (1984)
#' @name Pearce1984
#' @rdname Pearce1984
#' @description tectonic discrimination diagram for granites
#' @param Y vector with Y concentrations (ppm)
#' @param Nb vector with Nb concentrations (ppm)
#' @param Yb vector with Yb concentrations (ppm)
#' @param Ta vector with Ta concentrations (ppm)
#' @param Rb vector with Rb concentrations (ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param xlab (optional) x-axis label
#' @param ylab (optional) y-axis label
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @references Pearce, J. A., Harris, N. B. W., and Tindle, A. G.,
#'     1984, "Trace Element Discrimination Diagrams for the Tectonic
#'     Interpretation of Granitic Rocks": Journal of Petrology,
#'     v. 25, no. 4, p. 956-983.
#' @examples
#' data(test,package='GeoplotR')
#' @return a vector with tectonic affinities
NULL

#' @rdname Pearce1984
#' @examples
#' YNb(Y=test[,'Y'],Nb=test[,'Nb'])
#' @export
YNb <- function(Y=NULL,Nb=NULL,show.labels=TRUE,
                xlim=c(1,1000),ylim=c(1,1000),xlab,ylab,...){
    if (missing(xlab)) xlab <- 'Y (ppm)'
    if (missing(ylab)) ylab <- 'Nb (ppm)'
    invisible(xyplot(json=.YNb,X=Y,Y=Nb,show.labels=show.labels,
                     xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log='xy',...))
}

#' @rdname Pearce1984
#' @examples
#' YNbRb(Y=test[,'Y'],Nb=test[,'Nb'],Rb=test[,'Rb'])
#' @export
YNbRb <- function(Y=NULL,Nb=NULL,Rb=NULL,show.labels=FALSE,
                  xlim=c(1,1000),ylim=c(1,1000),xlab,ylab,...){
    if (missing(xlab)) xlab <- '(Y + Nb)(ppm)'
    if (missing(ylab)) ylab <- 'Rb (ppm)'
    invisible(xyplot(json=.YNbRb,X=Y+Nb,Y=Rb,show.labels=show.labels,
                     xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log='xy',...))    
}

#' @rdname Pearce1984
#' @examples
#' YbTa(Yb=test[,'Yb'],Ta=test[,'Ta'])
#' @export
YbTa <- function(Yb=NULL,Ta=NULL,show.labels=FALSE,
                 xlim=c(0.1,100),ylim=c(0.05,100),xlab,ylab,...){
    if (missing(xlab)) xlab <- 'Yb (ppm)'
    if (missing(ylab)) ylab <- 'Ta (ppm)'
    invisible(xyplot(json=.YbTa,X=Yb,Y=Ta,show.labels=show.labels,
                     xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log='xy',...))
}

#' @rdname Pearce1984
#' @examples
#' YbTaRb(Yb=test[,'Yb'],Ta=test[,'Ta'],Rb=test[,'Rb'])
#' @export
YbTaRb <- function(Yb=NULL,Ta=NULL,Rb=NULL,show.labels=TRUE,
                   xlim=c(0.9,110),ylim=c(1,1000),xlab,ylab,...){
    if (missing(xlab)) xlab <- '(Yb + Ta)(ppm)'
    if (missing(ylab)) ylab <- 'Rb (ppm)'
    invisible(xyplot(json=.YbTaRb,X=Yb+Ta,Y=Rb,show.labels=show.labels,
                     xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,log='xy',...))    
}
