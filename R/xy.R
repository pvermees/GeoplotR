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
#' LaYbTb(La_n=100,Yb_n=10)
#' @export
LaYb <- function(La_n,Yb_n,xlim=c(0,25),
                   ylim=c(0,160),show.labels=TRUE,...){
    xlab <- expression('Yb'[n])
    ylab <- expression('La'[n]*'/Yb'[n])
    invisible(xyplot(json=.LaYb,X=Yb_n,Y=La_n/Yb_n,xlim=xlim,ylim=ylim,
                     show.labels=show.labels,xlab=xlab,ylab=ylab,...))
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
CrY <- function(Cr,Y,xlim=c(1,100),ylim=c(1,10000),show.labels=TRUE,...){
    invisible(xyplot(json=.CrY,X=Cr,Y=Y,log='xy',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,xlab='Cr',ylab='Y',...))
}

#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with Na\eqn{_2}O concentrations (wt\%)
#' @param K2O vector with K\eqn{_2}O concentrations (wt\%)
#' @param SiO2 vector with SiO\eqn{_2} concentrations (wt\%)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param volcanic logical. Set to \code{FALSE} for plutonic rocks.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with rock types
#' @references Le Maitre, R. W., Streckeisen, A., Zanettin, B., Bas,
#'     M. J. L., Bonin, B., and Bateman, P., 2004,
#'     "Igneous Rocks: A Classification and Glossary of Terms":
#'     Cambridge University Press, v. -1, no. 70, p. 93–120.
#' @examples
#' data(test,package='GeoplotR')
#' TAS(test[,'Na2O'],test[,'K2O'],test[,'SiO2'])
#' @export
TAS <- function(Na2O=NULL,K2O=NULL,SiO2=NULL,
                xlim=c(35,80),ylim=c(0,15),
                volcanic=TRUE,short=TRUE,...){
    xlab <- expression('SiO'[2])
    ylab <- expression('Na'[2]*'O+K'[2]*'O')
    invisible(xyplot(json=.TAS,X=SiO2,Y=Na2O+K2O,
                     xlim=xlim,ylim=ylim,
                     short=short,xlab=xlab,ylab=ylab,...))
}

xyplot <- function(json,X=NULL,Y=NULL,xlim=range(X,na.rm=TRUE),
                   ylim=range(Y,na.rm=TRUE),xlab='x',ylab='y',
                   show.labels=FALSE,log='',short=FALSE,
                   test.polygons=FALSE,smooth=FALSE,...){
    graphics::plot(xlim,ylim,type='n',xlim=xlim,ylim=ylim,
                   xlab=xlab,ylab=ylab,bty='n',log=log)
    if (test.polygons){
        col <- 2
        for (pname in names(json$polygons)){
            xy <-  matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            graphics::polygon(xy,col=col)
            col <- col+1
        }
    } else {
        for (lname in names(json$lines)){
            xy <- matrix(unlist(json$lines[[lname]]),ncol=2,byrow=TRUE)
            if (smooth) shape <- 1
            else shape <- 0
            graphics::xspline(x=xy[,1],y=xy[,2],shape=shape,
                              lty=lty(json$line_type[[lname]]))
        }
    }
    if ( is.null(X) | is.null(Y) ){
        out <- NULL
    } else {
        pts <- cbind(X,Y)
        ns <- nrow(pts)
        out <- rep(NA,ns)
        col <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in 1:length(json$polygons)){
            pname <- pnames[i]
            pol <- matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            matched <- inside(pts=pts,pol=pol,log=log)
            out[matched] <- json$labels[[pname]]
            col[matched] <- i+1
        }
        graphics::points(pts,pch=21,bg=col,...)
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xy <- matrix(unlist(json$label_coords[[lname]]),ncol=2,byrow=TRUE)
            if (short) lab <- lname
            else lab <- json$labels[[lname]]
            a <- angle(json$angle[[lname]])
            graphics::text(xy,labels=lab,srt=a,pos=1)
        }
    }
    invisible(out)
}
