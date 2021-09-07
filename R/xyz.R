#' @title An-Ab-Or
#' @description An-Ab-Or diagram for tonalite, trondhjemite,
#'     granodiorite, and granite according to Barker (1979).
#' @param An vector with normative anorthite concentrations (wt\%)
#' @param Ab vector with normative albite concentrations (wt\%)
#' @param Or vector with normative orthoclase concentrations (wt\%)
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return a list with labels, indicating whether the input data are
#'     tonalite, trondhjemite, granodiorite, granite.
#' @references Barker, F., 1979, Trondhjemite: Definition,
#'     environment, and hypotheses of origin: p. 1-12, in Barker, F.,
#'     ed., Trondhjemites, Dacites, and Related Rocks, Elsevier,
#'     Amsterdam, 659 p.
#' @examples
#' AnAbOr(An=c(70,75,73),Ab=c(20,10,27),Or=c(10,15,0))
#' @export
AnAbOr <- function(An=NULL,Ab=NULL,Or=NULL,show.labels=TRUE,...){
    invisible(xyzplot(json=.AnAbOr,X=Ab,Y=An,Z=Or,
                      show.labels=show.labels,...))
}

#' @title Ti-Zr-Y
#' @description Ti-Zr-Y tectonic discrimination diagram
#' @param Ti vector with Ti concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param Y vector with Y concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Pearce'} for the nominal decision boundaries of Pearce
#'     and Cann (1973). The latter option has not been implemented
#'     yet.
#' @param plot either \code{'none'} to omit the plot, \code{'ternary'}
#'     for a ternary diagram, or \code{'logratio'} for a bivariate
#'     logratio plot
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
        out <- DA(uv=uv,da=da,ternary=ternary,f=c(1/100,1,3),pch=pch,bg=bg,...)
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

xyzplot <- function(json,X=NULL,Y=NULL,Z=NULL,f=rep(1,3),xyzlab=NULL,
                    show.labels=FALSE,short=FALSE,test.polygons=FALSE,
                    smooth=FALSE,pch=21,bg=NULL,dlwd=1,dcol='black',...){
    oldpar <- graphics::par(mar=rep(2,4),xpd=NA)
    if (is.null(xyzlab)) xyzlab <- json$axis
    ternaryplot(f=f,xyzlab=xyzlab,...)
    if (test.polygons){
        pcol <- 2
        for (pname in names(json$polygons)){
            xyz <-  matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            graphics::polygon(xyz2xy(xyz),col=pcol)
            pcol <- pcol+1
        }
    } else {
        for (lname in names(json$lines)){
            xyz <- matrix(unlist(json$lines[[lname]]),ncol=3,byrow=TRUE)
            xy <- xyz2xy(xyz)
            if (smooth) shape <- 1
            else shape <- 0
            graphics::xspline(x=xy[,1],y=xy[,2],shape=shape,border=dcol,
                              lty=lty(json$line_type[[lname]]),lwd=dlwd)
        }
    }
    if (is.null(X) | is.null(Y) | is.null(Z)){
        out <- NULL
    } else {
        XYZ <- sweep(cbind(X,Y,Z),2,f,'*')
        uv <- alr(XYZ,inverse=FALSE)
        XY <- xyz2xy(alr(uv,inverse=TRUE))
        ns <- nrow(XY)
        out <- rep(NA,ns)
        missingbg <- is.null(bg)
        if (missingbg) bg <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in 1:length(json$polygons)){
            pname <- pnames[i]
            xyz <- matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            matched <- inside(pts=XY,pol=xyz2xy(xyz))
            out[matched] <- ifelse(is.na(out[matched]),json$labels[[pname]],
                                   paste0(out[matched],' + ',json$labels[[pname]]))
            if (missingbg) bg[matched] <- i+1
        }
        ternarypoints(uv,pch=pch,bg=bg,...)
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xyz <- matrix(unlist(json$label_coords[[lname]]),ncol=3,byrow=TRUE)
            if (short) lab <- lname
            else lab <- json$labels[[lname]]
            graphics::text(xyz2xy(xyz),labels=lab,srt=json$angle[[lname]],pos=1)
        }
    }
    graphics::par(oldpar)
    invisible(out)
}
