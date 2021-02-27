#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with Na\eqn{_2}O concentrations (wt\%)
#' @param K2O vector with K\eqn{_2}O concentrations (wt\%)
#' @param SiO2 vector with SiO\eqn{_2} concentrations (wt\%)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param volcanic logical. Set to \code{FALSE} for plutonic rocks.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with rock types
#' @references Le Maitre, R. W., Streckeisen, A., Zanettin, B.,
#' Bas, M. J. L., Bonin, B., and Bateman, P., 2004,
#' "Igneous Rocks: A Classification and Glossary of Terms":
#' Cambridge University Press, v. -1, no. 70, p. 93–120.
#' @examples
#' data(test,package='GeoplotR')
#' TAS(test[,'Na2O'],test[,'K2O'],test[,'SiO2'])
#' @export
TAS <- function(Na2O=NULL,K2O=NULL,SiO2=NULL,xlim=c(35,80),
                ylim=c(0,15),volcanic=TRUE,...){
    xlab <- expression('SiO'[2])
    ylab <- expression('Na'[2]*'O+K'[2]*'O')
    invisible(xyplot(json=.TAS,X=SiO2,Y=Na2O+K2O,
                     xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...))
}

#' @title Y-Nb diagram
#' @description Y-Nb tectonic discrimination diagram for granites
#' @param Y vector with Y concentrations (ppm)
#' @param Nb vector with Nb concentrations (ppm)
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Pearce, J. A., Harris, N. B. W., and Tindle, A. G.,
#' 1984, "Trace Element Discrimination Diagrams
#' for the Tectonic Interpretation of Granitic Rocks":
#' Journal of Petrology, v. 25, no. 4, p. 956-983.
#' @examples
#' data(test,package='GeoplotR')
#' YNb(Y=test[,'Y'],Nb=test[,'Nb'])
#' @export
YNb <- function(Y=NULL,Nb=NULL,xlim=c(1,1000),ylim=c(1,1000),...){
    invisible(xyplot(json=.YNb,X=Y,Y=Nb,xlim=xlim,ylim=ylim,
                     xlab='Y (ppm)',ylab='Nb (ppm)',log='xy',...))
}
YNbRb <- function(Y=NULL,Nb=NULL,Rb=NULL,xlim=c(1,1000),ylim=c(1,1000),...){
    invisible(xyplot(json=.YNbRb,X=Y+Nb,Y=Rb,xlim=xlim,ylim=ylim,
                     xlab='(Y + Nb)(ppm)',ylab='Rb (ppm)',log='xy',...))    
}
YbTa <- function(Yb=NULL,Ta=NULL,xlim=c(0.1,100),ylim=c(0.1,100),...){
    invisible(xyplot(json=.YbTa,X=Yb,Y=Ta,xlim=xlim,ylim=ylim,
                     xlab='Yb (ppm)',ylab='Ta (ppm)',log='xy',...))
}
YbTaRb <- function(Yb=NULL,Ta=NULL,Rb=NULL,xlim=c(0.9,110),ylim=c(1,1000),...){
    invisible(xyplot(json=.YbTaRb,X=Yb+Ta,Y=Rb,xlim=xlim,ylim=ylim,
                     xlab='(Yb + Ta)(ppm)',ylab='Rb (ppm)',log='xy',...))    
}

xyplot <- function(json,X=NULL,Y=NULL,xlim=range(x,na.rm=TRUE),
                   ylim=range(y,na.rm=TRUE),xlab='x',ylab='y',
                   show.labels=FALSE,log='',...){
    graphics::plot(xlim,ylim,type='n',
                   xlab=xlab,ylab=ylab,bty='n',log=log)
    for (lname in names(json$lines)){
        xy <- matrix(unlist(json$lines[[lname]]),ncol=2,byrow=TRUE)
        graphics::lines(xy,lty=lty(json$line_type[[lname]]))
    }
    if ( is.null(X) | is.null(Y) ){
        out <- NULL
    } else {
        XY <- cbind(X,Y)
        ns <- nrow(XY)
        out <- rep(NA,ns)
        col <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in 1:length(json$polygons)){
            pname <- pnames[i]
            xy <- matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            matched <- apply(XY,MARGIN=1,FUN=inside,pol=xy)
            out[matched] <- json$labels[[pname]]
            col[matched] <- i+1
        }
        ternarypoints(XY,pch=21,bg=col,...)        
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xy <- matrix(unlist(json$label_coords[[lname]]),ncol=2,byrow=TRUE)
            graphics::text(xy,labels=json$labels[[lname]],
                           srt=json$angle[[lname]],pos=1)
        }
    }
    invisible(out)
}
