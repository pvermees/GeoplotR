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
TiV <- function(Ti=NULL,V=NULL,type=c('LDA','QDA','Shervais'),
                ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,short=TRUE,...){
    if (identical(type[1],'Shervais')){
        out <- TiV_nominal(Ti=Ti,V=V,pch=pch,bg=bg,
                           show.labels=show.labels,short=short,...)
    } else {
        uv <- alr(cbind(1e6-Ti-V,Ti,V))
        quadratic <- identical(type[1],'QDA')
        if (quadratic) da <- .TiV_QDA
        else da <- .TiV_LDA
        out <- DA(uv=uv,da=da,D2=TRUE,ternary=ternary,
                  xlab='Ti',ylab='V',f=c(1,100,5000),pch=pch,bg=bg,...)
        plotlabels(diagram='TiV',ternary=ternary,f=c(1,1,5000),
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
}
TiV_nominal <- function(Ti=NULL,V=NULL,pch=21,bg=NULL,
                        show.labels=TRUE,short=TRUE,...){
    invisible(xyplot(json=.TiV_nominal,X=Ti/1000,Y=V,pch=pch,bg=bg,
                     short=short,show.labels=show.labels,...))
}
