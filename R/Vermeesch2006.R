#' @title Ti-Zr-Y
#' @description Ti-Zr-Y tectonic discrimination diagram
#' @param Ti vector with Ti concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param Y vector with Y concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Pearce'} for the nominal decision boundaries of Pearce
#'     and Cann (1973).
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with labels for \code{MORB}, \code{IAB} and
#'     \code{OIB}.
#' @references Pearce, J. A., and Cann, J. R., 1973, Tectonic setting
#'     of basic volcanic rocks determined using trace element
#'     analyses: Earth and Planetary Science Letters, v. 19, no. 2,
#'     p. 290-300.
#' @examples
#' data(test,package='GeoplotR')
#' TiZrY(Ti=wtpct2ppm(test[,'TiO2']),
#'       Zr=test[,'Zr'],Y=test[,'Y'],
#'       type='QDA',plot='ternary')
#' @export
TiZrY <- function(Ti=NULL,Zr=NULL,Y=NULL,
                  type=c('LDA','QDA','Pearce'),
                  ternary=TRUE,pch=21,bg=NULL,
                  show.labels=FALSE,short=TRUE,...){
    if (identical(type[1],'Pearce')) {
        out <- TiZrY_nominal(Ti=Ti,Zr=Zr,Y=Y,pch=pch,bg=bg,
                             show.labels=show.labels,short=short,...)
    } else {
        uv <- alr(cbind(Ti,Zr,Y))
        quadratic <- identical(type[1],'QDA')
        if (quadratic) da <- .TiZrY_QDA
        else da <- .TiZrY_LDA
        out <- DA(uv=uv,da=da,ternary=ternary,f=c(1/100,1,3),
                  xlab='Ti',ylab='Zr',zlab='Y',pch=pch,bg=bg,...)
        plotlabels(diagram='TiZrY',ternary=ternary,f=c(1/100,1,3),
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
    invisible(out)
}
TiZrY_nominal <- function(Ti=NULL,Zr=NULL,Y=NULL,pch=21,bg=NULL,
                          show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.TiZrY_nominal,X=Ti,Y=Zr,Z=Y,f=c(0.01,1,3),
                      pch=pch,bg=bg,short=short,show.labels=show.labels,...))
}

#' @title Ti-V
#' @description Ti-V tectonic discrimination diagram
#' @param Ti vector with Ti concentrations (ppm)
#' @param V vector with V concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Shervais'} for the nominal decision boundaries of
#'     Shervais (1982).
#' @param ternary logical. If \code{TRUE}, produces a ternary diagram.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with labels for \code{MORB}, \code{IAB} and
#'     \code{OIB}.
#' @references Shervais, J.W., 1982. Ti-V plots and the petrogenesis
#'     of modern and ophiolitic lavas. Earth and Planetary Science
#'     Letters, 59(1), pp.101-118.
#' @examples
#' data(test,package='GeoplotR')
#' TiV(Ti=wtpct2ppm(test[,'TiO2']),V=test[,'V'],type='Shervais')
#' @export
TiV <- function(Ti=NULL,V=NULL,type=c('Shervais','LDA','QDA'),
                ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,
                short=TRUE,xlim=NULL,ylim=NULL,...){
    if (identical(type[1],'Shervais')){
        out <- TiV_nominal(Ti=Ti,V=V,pch=pch,bg=bg,show.labels=show.labels,
                           short=short,xlim=xlim,ylim=ylim,...)
    } else {
        tot <- 1e6
        uv <- alr(cbind(tot-Ti-V,Ti,V))
        quadratic <- identical(type[1],'QDA')
        if (quadratic) da <- .TiV_QDA
        else da <- .TiV_LDA
        out <- DA(uv=uv,da=da,D2=TRUE,tot=tot,ternary=ternary,
                  xlab='Ti',ylab='V',f=c(1,100,5000),
                  pch=pch,bg=bg,xlim=xlim,ylim=ylim,...)
        plotlabels(diagram='TiV',ternary=ternary,f=c(1,100,5000),linear=TRUE,
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
    invisible(out)
}
TiV_nominal <- function(Ti=NULL,V=NULL,pch=21,bg=NULL,show.labels=TRUE,
                        short=TRUE,xlim=NULL,ylim=NULL,...){
    good <- !(is.na(Ti) | is.na(V))
    if (is.null(xlim)) xlim <- getlimits(x=Ti[good],m=0,M=25)/1000
    if (is.null(ylim)) ylim <- getlimits(x=V[good],m=0,M=600)
    invisible(xyplot(json=.TiV_nominal,X=Ti/1000,Y=V,pch=pch,bg=bg,
                     short=short,show.labels=show.labels,
                     xlim=xlim,ylim=ylim,...))
}

#' @title Ti-Zr
#' @description Ti-Zr tectonic discrimination diagram.
#' @param Ti vector with Ti concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis,
#'     \code{'Pearce'} for the nominal decision boundaries of Pearce
#'     and Cann (1973), or \code{'Dilek'} for the nominal decision
#'     boundaries of Dilek and Furnes (2009).
#' @param ternary logical. If \code{TRUE}, produces a ternary diagram.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Pearce, J.A. and Cann, J.R., 1973. Tectonic setting of
#'     basic volcanic rocks determined using trace element
#'     analyses. Earth and Planetary Science Letters, 19(2),
#'     pp. 290-300.
#' 
#' Dilek, Y. and Furnes, H., 2009, Structure and geochemistry of
#' Tethyan ophiolites and their petrogenesis in subduction rollback
#' systems.  Lithos, v. 113, p. 1-20.
#' @examples
#' data(test,package='GeoplotR')
#' ZrTi(Zr=c(50,40,60),Ti=c(1000,4000,6000))
#' @export
ZrTi <- function(Zr=NULL,Ti=NULL,type=c('LDA','QDA','Pearce','Dilek'),
                 ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,
                 short=TRUE,xlim=NULL,ylim=NULL,...){
    if (identical(type[1],'Pearce')){
        out <- ZrTi_nominal(Zr=Zr,Ti=Ti,pch=pch,bg=bg,show.labels=show.labels,
                            short=short,xlim=xlim,ylim=ylim,...)
    } else if (identical(type[1],'Dilek')){
        out <- TiZr_Dilek(Ti=Ti,Zr=Zr,pch=pch,bg=bg,show.labels=show.labels,
                          short=short,xlim=xlim,ylim=ylim,...)        
    } else {
        tot <- 1e6
        uv <- alr(cbind(tot-Ti-Zr,Zr,Ti))
        quadratic <- identical(type[1],'QDA')
        if (quadratic) da <- .ZrTi_QDA
        else da <- .ZrTi_LDA
        out <- DA(uv=uv,da=da,D2=TRUE,tot=tot,ternary=ternary,
                  xlab='Zr',ylab='Ti',f=c(1,15000,200),
                  pch=pch,bg=bg,xlim=xlim,ylim=ylim,...)
        plotlabels(diagram='ZrTi',ternary=ternary,f=c(1,15000,200),linear=TRUE,
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
    invisible(out)
}
ZrTi_nominal <- function(Zr=NULL,Ti=NULL,xlim=NULL,ylim=NULL,
                         show.labels=TRUE,short=FALSE,...){
    good <- !(is.na(Zr) | is.na(Ti))
    if (is.null(xlim)) xlim <- getlimits(x=Zr[good],m=0,M=200)
    if (is.null(ylim)) ylim <- getlimits(x=Ti[good],m=0,M=15000)
    invisible(xyplot(json=.ZrTi_nominal,X=Zr,Y=Ti,log='',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,short=short,smooth=FALSE,...))
}
TiZr_Dilek <- function(Ti=NULL,Zr=NULL,xlim=NULL,ylim=NULL,
                       show.labels=TRUE,short=FALSE,...){
    good <- !(is.na(Ti) | is.na(Zr))
    if (is.null(xlim)) xlim <- getlimits(x=Ti[good],m=0,M=10000)
    if (is.null(ylim)) ylim <- getlimits(x=Zr[good],m=0,M=150)
    invisible(xyplot(json=.TiZr,X=Ti,Y=Zr,log='',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,short=short,smooth=TRUE,...))
}
