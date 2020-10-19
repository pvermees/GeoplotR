#' @title TAS diagram
#' @description Total Alkali Silica diagram
#' @param Na2O vector with \eqn{Na_2O} concentrations (wt\%)
#' @param K2O vector with \eqn{K_2O} concentrations (wt\%)
#' @param SiO2 vector with \eqn{SiO_2} concentrations (wt\%)
#' @param plot logical. If \code{FALSE}, suppresses the graphical
#'     output.
#' @param labels logical. Label to mark the different fields of the
#'     TAS diagram.
#' @param ... additional arguments to the generic \code{points}
#'     function.
#' @return a vector with rock types
#' @examples
#' data(test,package='GeoplotR')
#' TAS(Na2O=test[,'NA2O(WT%)'],
#'     K2O=test[,'K2O(WT%)'],
#'     SiO2=test[,'SIO2(WT%)'])
#' @export
TAS <- function(Na2O,K2O,SiO2,plot=TRUE,labels=FALSE,...){
    plot(c(35,90),c(0,20),type='n')
    nc <- length(tas$cords)
    labs <- names(tas$cords)
    TA <- Na2O+K2O
    S <- SiO2
    graphics::points(S,TA,...)
    ns <- length(S)
    out <- rep(NA,ns)
    for (i in 1:nc){
        cord <- unlist(tas$cords[i])
        xy <- matrix(cord,ncol=2,byrow=TRUE)
        graphics::lines(xy)
        if (labels){
            xyl <- colMeans(xy)
            graphics::text(x=xyl[1],y=xyl[2],labels=labs[i])
        }
        for (j in 1:ns){
            good <- !(is.na(S[j]) | is.na(TA[j]))
            if (good && inside(x=S[j],y=TA[j],X=xy[,1],Y=xy[,2])){
                out[j] <- tas$Plutonic[[labs[i]]]
            }
        }
    }
    invisible(out)
}

# x,y: coordinates of the sample
# X,Y: coordinates of the class
inside <- function(x,y,X,Y){
    ch <- grDevices::chull(x=X,y=Y)
    CH <- grDevices::chull(x=c(X,x),y=c(Y,y))
    if (identical(ch,CH)) return(TRUE)
    else return(FALSE)
}
