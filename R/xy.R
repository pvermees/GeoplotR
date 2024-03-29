#' @title La/Yb-Yb
#' @description (La/Yb)_N vs. Yb_N diagram of Martin (1986), which
#'     discriminates between Archean TTG suite and more modern
#'     adakites on the one hand, and classical island arcs on the
#'     other hand.
#' @param La_n vector with chondritic normalised La concentrations
#' @param Yb_n vector with chondritic normalised Yb concentrations
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Martin, H., 1986. Effect of steeper Archean geothermal
#'     gradient on geochemistry of subduction-zone magmas. Geology,
#'     14(9), pp.753-756.
#' @examples
#' data(test,package='GeoplotR')
#' LaYb(La_n=100,Yb_n=10)
#' @export
LaYb <- function(La_n=NULL,Yb_n=NULL,xlim=NULL,ylim=NULL,show.labels=TRUE,...){
    good <- !(is.na(La_n) | is.na(Yb_n))
    if (is.null(xlim)) xlim <- getlimits(x=Yb_n[good],m=0,M=20)
    if (is.null(ylim)) ylim <- getlimits(x=La_n[good]/Yb_n[good],m=0,M=160)
    xlab <- expression('Yb'[n])
    ylab <- expression('La'[n]*'/Yb'[n])
    invisible(xyplot(json=.LaYb,X=Yb_n,Y=La_n/Yb_n,xlim=xlim,ylim=ylim,
                     show.labels=show.labels,xlab=xlab,ylab=ylab,...))
}

#' @title Sr/Y-Y
#' @description (Sr/Y) vs. Y diagram of Drummond (1990), which
#'     discriminates between Archean TTG suite and more modern
#'     adakites on the one hand, and classical island arcs on the
#'     other hand.
#' @param Sr vector with Sr concentrations (in ppm)
#' @param Y vector with Y concentrations (in ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Drummond MS, Defant MJ (1990) A model for
#'     trondhjemite-tonalite-dacite genesis and crustal growth via
#'     slab melting: Archean to modern comparisons. J Geophys Res
#'     95:21503– 21521
#' @examples
#' data(test,package='GeoplotR')
#' SrY(Sr=1000,Y=10)
#' @export
SrY <- function(Sr=NULL,Y=NULL,xlim=NULL,ylim=NULL,show.labels=TRUE,...){
    good <- !(is.na(Sr) | is.na(Y))
    if (is.null(xlim)) xlim <- getlimits(x=Y[good],m=0,M=60)
    if (is.null(ylim)) ylim <- getlimits(x=Sr[good]/Y[good],m=0,M=500)
    invisible(xyplot(json=.SrY,X=Y,Y=Sr/Y,xlim=xlim,ylim=ylim,
                     show.labels=show.labels,...))
}

#' @title Cr-Y
#' @description Cr vs Y discriminant for MORB, IAT, and boninites of
#'     Dilek et al. (2007)
#' @param Cr vector with Cr concentrations (ppm)
#' @param Y vector with Y concentrations (ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Dilek, Y., Furnes, H., and Shallo, M., 2007,
#'     Suprasubduction zone ophiolite formation along the periphery of
#'     Mesozoic Gondwana. Gondwana Research, v. 11, p. 453-475.
#' @examples
#' data(test,package='GeoplotR')
#' CrY(Cr=test[,'Cr'],Y=test[,'Y'])
#' @export
CrY <- function(Cr=NULL,Y=NULL,xlim=NULL,ylim=NULL,show.labels=TRUE,...){
    good <- !(is.na(Cr) | is.na(Y))
    if (is.null(xlim)) xlim <- getlimits(x=Cr[good],m=1,M=100)
    if (is.null(ylim)) ylim <- getlimits(x=Y[good],m=1,M=10000)
    invisible(xyplot(json=.CrY,X=Cr,Y=Y,log='xy',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,smooth=TRUE,...))
}

#' @title Th-Co
#' @description Th vs Co discrimination diagram of Hastie et
#'     al. (2007) for tholeiitic, calc-alkaline, and more alkaline arc
#'     series, and basaltic, andesitic, and rhyolitic compositions.
#' @param Th vector with Th concentrations (ppm)
#' @param Co vector with Co concentrations (ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Hastie, A.R., Kerr, A.C., Pearce, J.A., Mitchell, S.F.,
#'     2007, Classification of Altered Volcanic Island Arc Rocks using
#'     Immobile Trace Elements: Development of the Th-Co
#'     Discrimination Diagram.  Journal of Petrology, v. 48,
#'     p. 2341-2357.
#' @examples
#' data(test,package='GeoplotR')
#' ThCo(Th=test[,'Th'],Co=test[,'Co'])
#' @export
ThCo <- function(Th=NULL,Co=NULL,xlim=NULL,ylim=NULL,show.labels=TRUE,...){
    good <- !(is.na(Th) | is.na(Co))
    if (is.null(xlim)) xlim <- getlimits(x=Co[good],m=70,M=0)
    if (is.null(ylim)) ylim <- getlimits(x=Th[good],m=0.01,M=100)
    invisible(xyplot(json=.ThCo,X=Co,Y=Th,log='y',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,smooth=TRUE,...))
}

#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with Na\eqn{_2}O concentrations (wt\%)
#' @param K2O vector with K\eqn{_2}O concentrations (wt\%)
#' @param SiO2 vector with SiO\eqn{_2} concentrations (wt\%)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param plutonic logical. If \code{TRUE}, uses the nomenclature of
#'     Middlemost (1994).
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with rock types
#' @references Le Maitre, R. W., Streckeisen, A., Zanettin, B., Bas,
#'     M. J. L., Bonin, B., and Bateman, P., 2004,
#'     "Igneous Rocks: A Classification and Glossary of Terms":
#'     Cambridge University Press, v. -1, no. 70, p. 93–120.
#'
#' Middlemost, E.A., 1994. Naming materials in the magma/igneous rock
#' system. Earth-science reviews, 37(3-4), pp.215-224.
#' @examples
#' data(test,package='GeoplotR')
#' TAS(test[,'Na2O'],test[,'K2O'],test[,'SiO2'])
#' @export
TAS <- function(Na2O=NULL,K2O=NULL,SiO2=NULL,
                xlim=NULL,ylim=NULL,plutonic=FALSE,
                show.labels=TRUE,short=TRUE,...){
    good <- !(is.na(Na2O) | is.na(K2O) | is.na(SiO2))
    if (is.null(xlim)) xlim <- getlimits(x=SiO2[good],m=35,M=80)
    if (is.null(ylim)) ylim <- getlimits(x=Na2O[good]+K2O[good],m=0,M=15)
    xlab <- expression('SiO'[2])
    ylab <- expression('Na'[2]*'O+K'[2]*'O')
    if (plutonic) json <- .TAS_plutonic
    else json <- .TAS
    invisible(xyplot(json=json,X=SiO2,Y=Na2O+K2O,
                     xlim=xlim,ylim=ylim,show.labels=show.labels,
                     short=short,xlab=xlab,ylab=ylab,...))
}

xyplot <- function(json,X=NULL,Y=NULL,xlim=range(X,na.rm=TRUE),
                   ylim=range(Y,na.rm=TRUE),xlab=NULL,ylab=NULL,
                   show.labels=FALSE,log='',short=FALSE,
                   test.polygons=FALSE,smooth=FALSE,bg=NULL,
                   pch=21,dlwd=1,dcol='black',...){
    if (is.null(xlab)) xlab <- json$axis[1]
    if (is.null(ylab)) ylab <- json$axis[2]
    graphics::plot(xlim,ylim,type='n',xlim=xlim,ylim=ylim,
                   xlab=xlab,ylab=ylab,bty='n',log=log)
    if (test.polygons){
        pcol <- 2
        for (pname in names(json$polygons)){
            xy <-  matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            graphics::polygon(xy,col=pcol)
            pcol <- pcol+1
        }
    } else {
        for (lname in names(json$lines)){
            xy <- matrix(unlist(json$lines[[lname]]),ncol=2,byrow=TRUE)
            if (smooth) shape <- 1
            else shape <- 0
            graphics::xspline(x=xy[,1],y=xy[,2],shape=shape,col=dcol,
                              lwd=dlwd,lty=lty(json$line_type[[lname]]))
        }
    }
    if ( is.null(X) | is.null(Y) ){
        out <- NULL
    } else {
        pts <- cbind(X,Y)
        ns <- nrow(pts)
        out <- rep(NA,ns)
        missingbg <- is.null(bg)
        if (missingbg) bg <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in 1:length(json$polygons)){
            pname <- pnames[i]
            pol <- matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            matched <- inside(pts=pts,pol=pol,log=log)
            out[matched] <- ifelse(is.na(out[matched]),json$labels[[pname]],
                                   paste0(out[matched],' + ',json$labels[[pname]]))
            if (missingbg) bg[matched] <- i+1
        }
        graphics::points(pts,bg=bg,pch=pch,...)
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xy <- matrix(unlist(json$label_coords[[lname]]),ncol=2,byrow=TRUE)
            if (short) lab <- lname
            else lab <- json$labels[[lname]]
            a <- angle(json$angle[[lname]])
            graphics::text(xy,labels=lab,srt=a,xpd=TRUE)
        }
    }
    invisible(out)
}
