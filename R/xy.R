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
YNb <- function(Y=NULL,Nb=NULL,xlim=c(1,1000),ylim=c(1,1000),...){
    invisible(xyplot(json=.YNb,X=Y,Y=Nb,xlim=xlim,ylim=ylim,
                     xlab='Y (ppm)',ylab='Nb (ppm)',log='xy',...))
}

#' @rdname Pearce1984
#' @examples
#' YNbRb(Y=test[,'Y'],Nb=test[,'Nb'],Rb=test[,'Rb'])
#' @export
YNbRb <- function(Y=NULL,Nb=NULL,Rb=NULL,xlim=c(1,1000),ylim=c(1,1000),...){
    invisible(xyplot(json=.YNbRb,X=Y+Nb,Y=Rb,xlim=xlim,ylim=ylim,
                     xlab='(Y + Nb)(ppm)',ylab='Rb (ppm)',log='xy',...))    
}

#' @rdname Pearce1984
#' @examples
#' YbTa(Yb=test[,'Yb'],Ta=test[,'Ta'])
#' @export
YbTa <- function(Yb=NULL,Ta=NULL,xlim=c(0.1,100),ylim=c(0.05,100),...){
    invisible(xyplot(json=.YbTa,X=Yb,Y=Ta,xlim=xlim,ylim=ylim,
                     xlab='Yb (ppm)',ylab='Ta (ppm)',log='xy',...))
}

#' @rdname Pearce1984
#' @examples
#' YbTaRb(Yb=test[,'Yb'],Ta=test[,'Ta'],Rb=test[,'Rb'])
#' @export
YbTaRb <- function(Yb=NULL,Ta=NULL,Rb=NULL,xlim=c(0.9,110),ylim=c(1,1000),...){
    invisible(xyplot(json=.YbTaRb,X=Yb+Ta,Y=Rb,xlim=xlim,ylim=ylim,
                     xlab='(Yb + Ta)(ppm)',ylab='Rb (ppm)',log='xy',...))    
}

xyplot <- function(json,X=NULL,Y=NULL,xlim=range(X,na.rm=TRUE),
                   ylim=range(Y,na.rm=TRUE),xlab='x',ylab='y',
                   show.labels=FALSE,log='',short=FALSE,
                   test.polygons=FALSE,smooth=FALSE,...){
    graphics::plot(xlim,ylim,type='n',
                   xlab=xlab,ylab=ylab,bty='n',log=log)
    if (test.polygons){
        for (pname in names(json$polygons)){
            xy <-  matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            graphics::polygon(xy)
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
            matched <- inside(pts=pts,pol=pol)
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
            graphics::text(xy,labels=lab,srt=json$angle[[lname]],pos=1)
        }
    }
    invisible(out)
}
