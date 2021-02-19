#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with Na\eqn{_2}O concentrations (wt\%)
#' @param K2O vector with K\eqn{_2}O concentrations (wt\%)
#' @param SiO2 vector with SiO\eqn{_2} concentrations (wt\%)
#' @param plot logical. If \code{FALSE}, suppresses the graphical
#'     output.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param labels logical. Label to mark the different fields of the
#'     TAS diagram.
#' @param volcanic logical. Set to \code{FALSE} for plutonic rocks.
#' @param ... additional arguments to the generic \code{points}
#'     function.
#' @return a vector with rock types
#' @examples
#' data(test,package='GeoplotR')
#' TAS(test[,'Na2O'],test[,'K2O'],test[,'SiO2'])
#' @export
TAS <- function(Na2O,K2O,SiO2,plot=TRUE,xlim=c(35,90),
                ylim=c(0,20),labels=FALSE,volcanic=TRUE,...){
    if (plot){
        xlab <- expression('SiO'[2])
        ylab <- expression('Na'[2]*'O+K'[2]*'O')
        if (labels) tags <- names(.TAS$cords)
        else tags <- NULL
        xyplot(x=SiO2,y=Na2O+K2O,xlim=xlim,ylim=ylim,
               xlab=xlab,ylab=ylab,coords=.TAS$coords,
               tags=tags,...)
    }
    if (volcanic) classes <- .TAS$Volcanic
    else classes <- .TAS$Plutonic
    out <- xyclassify(x=SiO2,y=Na2O+K2O,coords=.TAS$coords,classes=classes)
    invisible(out)
}

PHT84 <- function(Nb,Y,Yb,Ta,Rb,plot=TRUE,labels=FALSE,...){
    out <- if (missing(Rb) & missing(Ta) & missing(Yb))
        NbY(Nb=Nb,Y=Y,plot=plot,labels=labels,...)
    else if (missing(Rb) & missing(Y) & missing(Nb))
        TaYb(Ta=Ta,Yb=Yb,plot=plot,labels=labels,...)
    else if (missing(Yb) & missing(Ta))
        RbYNb(Rb=Rb,Y=Y,Nb=Nb,plot=plot,labels=labels,...)
    else if (missing(Rb) | missing(Yb) | missing(Ta))
        stop('You must supply either Nb & Y; Ta & Nb; Rb, Y & Nb or Rb, Yb & Ta.')
    else
        RbYbTa(Rb=Rb,Yb=Yb,Ta=Ta,plot=plot,labels=labels,...)
    invisible(out)
}
NbY <- function(Nb,Y,plot=TRUE,labels=FALSE,...){
    invisible(xyhelper(x=Nb,y=Y,coords=.PHT84$coords2,
                       plot=plot,xlim=c(1,1000),ylim=c(1,1000),
                       labels=labels,log='xy',...))
}
TaYb <- function(Ta,Yb,plot=TRUE,labels=FALSE,...){
    invisible(xyhelper(x=Ta,y=Yb,coords=.PHT84$coords3,
                       plot=plot,xlim=c(0.1,100),ylim=c(0.1,100),
                       labels=labels,log='xy',...))
}
RbYNb <- function(Rb,Y,Nb,plot=TRUE,labels=FALSE,...){
    invisible(xyhelper(x=Rb,y=Y+Nb,coords=.PHT84$coords0,
                       plot=plot,xlim=c(1,1000),ylim=c(1,1000),
                       labels=labels,log='xy',...))
}
RbYbTa <- function(Rb,Yb,Ta,plot=TRUE,labels=FALSE,...){
    invisible(xyhelper(x=Rb,y=Yb+Ta,coords=.PHT84$coords1,
                       plot=plot,xlim=c(0.5,200),ylim=c(1,1000),
                       labels=labels,log='xy',...))
}

xyhelper <- function(x,y,coords,plot=TRUE,labels=FALSE,log='',...){
    if (plot){
        if (labels) tags <- coords$Labels
        else tags <- NULL
        xyplot(x=x,y=y,xlab=coords$xLabel,ylab=coords$yLabel,
               coords=coords$BaseLines,tags=tags,log=log,...)
    }
    classes <- coords$Labels
    out <- xyclassify(x=x,y=y,coords=coords$BaseLines,classes=classes)
    invisible(out)
}

xyplot <- function(x,y,xlim=range(x,na.rm=TRUE),ylim=range(y,na.rm=TRUE),
                   xlab='x',ylab='y',coords,tags=NULL,log='',...){
    graphics::plot(xlim,ylim,type='n',xlab=xlab,ylab=ylab,bty='n',log=log)
    for (i in 1:length(coords)){
        coord <- unlist(coords[i])
        xy <- matrix(coord,ncol=2,byrow=TRUE)
        graphics::lines(xy,col='gray50')
        if (!is.null(tags)){
            xyl <- colMeans(xy)
            graphics::text(x=xyl[1],y=xyl[2],labels=tags[i],col='gray50')
        }
    }
    graphics::points(x,y,...)
}

xyclassify <- function(x,y,coords,classes){
    nc <- length(coords)
    ns <- length(x)
    out <- rep(NA,ns)
    tags <- names(coords)
    for (i in 1:nc){
        coord <- unlist(coords[i])
        xy <- matrix(coord,ncol=2,byrow=TRUE)
        for (j in 1:ns){
            if (!(is.na(x[j]) | is.na(y[j]))){
                if (inside(x=x[j],y=y[j],X=xy[,1],Y=xy[,2])){
                    out[j] <- classes[[tags[i]]]
                }
            }
        }
    }
    out
}

# x,y: coordinates of the sample
# X,Y: coordinates of the class
inside <- function(x,y,X,Y){
    ch <- grDevices::chull(x=X,y=Y)
    CH <- grDevices::chull(x=c(X,x),y=c(Y,y))
    if (identical(ch,CH)) return(TRUE)
    else return(FALSE)
}
