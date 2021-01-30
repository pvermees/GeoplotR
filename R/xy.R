#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with Na\eqn{_2}O concentrations (wt\%)
#' @param K2O vector with K\eqn{_2}O concentrations (wt\%)
#' @param SiO2 vector with SiO\eqn{_2} concentrations (wt\%)
#' @param plot logical. If \code{FALSE}, suppresses the graphical
#'     output.
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
TAS <- function(Na2O,K2O,SiO2,plot=TRUE,labels=FALSE,volcanic=TRUE,...){
    if (plot){
        xlim <- c(35,90)
        ylim <- c(0,20)
        xlab <- expression('SiO'[2])
        ylab <- expression('Na'[2]*'O+K'[2]*'O')
        if (labels) tags <- names(.TAS$cords)
        else tags <- NULL
        xyplot(x=SiO2,y=Na2O+K2O,xlim=xlim,ylim=ylim,
               xlab=xlab,ylab=ylab,cords=.TAS$cords,
               tags=tags,...)
    }
    if (volcanic) classes <- .TAS$Volcanic
    else classes <- .TAS$Plutonic
    out <- xyclassify(x=SiO2,y=Na2O+K2O,cords=.TAS$cords,classes=classes)
    invisible(out)
}

xyplot <- function(x,y,xlim=range(x),ylim=range(y),
                   xlab='x',ylab='y',cords,tags=NULL,...){
    graphics::plot(xlim,ylim,type='n',xlab=xlab,ylab=ylab,bty='n')
    for (i in 1:length(cords)){
        cord <- unlist(cords[i])
        xy <- matrix(cord,ncol=2,byrow=TRUE)
        graphics::lines(xy,col='gray50')
        if (!is.null(tags)){
            xyl <- colMeans(xy)
            graphics::text(x=xyl[1],y=xyl[2],labels=tags[i],col='gray50')
        }
    }
    graphics::points(x,y,...)
}

xyclassify <- function(x,y,cords,classes){
    nc <- length(cords)
    ns <- length(x)
    out <- rep(NA,ns)
    tags <- names(cords)
    for (i in 1:nc){
        cord <- unlist(cords[i])
        xy <- matrix(cord,ncol=2,byrow=TRUE)
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
