# aa*x + bb*y + cc = 0
point2line <- function(xy,aa,bb,cc){
    x <- xy[1]
    y <- xy[2]
    x0 <- (bb*(bb*x-aa*y)-aa*cc)/(aa^2+bb^2)
    y0 <- (aa*(aa*y-bb*x)-bb*cc)/(aa^2+bb^2)
    c(x0,y0)
}

nearestpoint <- function(xy,p){
    xy1 <- point2line(xy,p[3],-1,p[2]-p[1]*p[3])
    xy2 <- point2line(xy,p[4],-1,p[2]-p[1]*p[4])
    offside1 <- (xy1[2]<p[2])
    offside2 <- (xy2[1]<p[1])
    if (offside1){
        if (offside2) out <- p[1:2]
        else out <- xy2
    } else if (offside2){
        out <- xy1
    } else {
        d1 <- sqrt(sum((xy-xy1)^2))
        d2 <- sqrt(sum((xy-xy2)^2))
        if (d1<d2) out <- xy1
        else out <- xy2
    }
    out
}

#' @title Bowen-Fenner index
#' @description Bowen-Fenner index for tholeiitic and calc-alkaline
#'     magma series
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param F vector with (FeO + FeO\eqn{_2}O\eqn{_3}) concentrations
#'     (in wt\%)
#' @param M vector with MgO concentrations (in wt\%)
#' @param T vector with TiO\eqn{_2} concentrations (in wt\%)
#' @param twostage logical. If \code{TRUE}, uses the two-stage magma
#'     evolution model. If \code{TRUE}, uses a simple one-stage
#'     model. Only used if \code{T} has been omitted.
#' @return the Bowen-Fenner indices of the samples, where positive
#'     values indicate calc-alkaline and negative numbers tholeiitic
#'     compositions.
#'
#' If \code{F} is provided, \code{T} is omitted, and \code{twostage}
#' is \code{FALSE}: returns the \eqn{BF_1} index of Vermeesch and
#' Pease (2021).
#'
#' If \code{F} is provided, \code{T} is omitted, and \code{twostage}
#' is \code{TRUE}: returns the \eqn{BF_2^F} index.
#'
#' If \code{F} is omitted and \code{T} is provided, returns the
#' \eqn{BF_2^T} index.
#'
#' If both \code{F} and \code{T} are provided, returns the \eqn{BF_2}
#' index, which is the mean of \eqn{BF_2^F} and \eqn{BF_2^T}.
#'
#' @references Vermeesch, P. and Pease, V. 2021, A genetic
#'     classification of the tholeiitic and calc-alkaline magma
#'     series, Geochemical Perspective Letters.
#' @examples
#' data(cath,package='GeoplotR')
#' 
#' oldpar <- par(mfrow=c(2,2),mar=c(4,4.5,2,0))
#' on.exit(par(oldpar))
#' 
#' A <- cath[,'Na2O']+cath[,'K2O']
#' F <- cath[,'FeOT']
#' T <- cath[,'TiO2']
#' M <- cath[,'MgO']
#' 
#' CA <- (cath$affinity=='ca')
#' TH <- (cath$affinity=='th')
#' 
#' dCA <- list()
#' dTH <- list()
#' 
#' titles <- c('BF1','BF2F','BF2T','BF2')
#' 
#' dCA[[1]] <- density(BF(A=A[CA],F=F[CA],M=M[CA],twostage=FALSE))
#' dTH[[1]] <- density(BF(A=A[TH],F=F[TH],M=M[TH],twostage=FALSE))
#' 
#' dCA[[2]] <- density(BF(A=A[CA],F=F[CA],M=M[CA],twostage=TRUE))
#' dTH[[2]] <- density(BF(A=A[TH],F=F[TH],M=M[TH],twostage=TRUE))
#' 
#' dCA[[3]] <- density(BF(A=A[CA],M=M[CA],T=T[CA],twostage=TRUE))
#' dTH[[3]] <- density(BF(A=A[TH],M=M[TH],T=T[TH],twostage=TRUE))
#' 
#' dCA[[4]] <- density(BF(A=A[CA],F=F[CA],M=M[CA],T=T[CA],twostage=TRUE))
#' dTH[[4]] <- density(BF(A=A[TH],F=F[TH],M=M[TH],T=T[TH],twostage=TRUE))
#' 
#' for (i in 1:4){
#'     matplot(cbind(dCA[[i]]$x,dTH[[i]]$x),
#'             cbind(dCA[[i]]$y,dTH[[i]]$y),
#'             type='l',lty=1,col=c('red','black'),
#'             xlab='BF',ylab='')
#'     title(titles[i])
#'     legend('topleft',legend=c('Cascades','Iceland'),
#'            lty=1,col=c('red','black'))
#' }
#' 
#' @export
BF <- function(A,F,M,T,twostage=TRUE){
    if (missing(F) & missing(T)){
        stop('The Bowen-Fenner index requires either F or T.')
    } else if (missing(T)){
        if (twostage){
            out <- BF_Fe2(A,F,M)
        } else {
            out <- BF_Fe1(A,F,M)
        }
    } else if (missing(F)){
        out <- BF_Ti(A,T,M)
    } else {
        out <- (BF_Fe2(A,F,M) + BF_Ti(A,T,M))/2
    }
    out
}
BF_Fe1 <- function(A,F,M,project=FALSE){
    uv <- alr(cbind(F,A,M))
    BF_helper(uv,pars=.BF_Fe1,twostage=FALSE)
}
BF_Fe2 <- function(A,F,M,project=FALSE){    
    uv <- alr(cbind(F,A,M))
    BF_helper(uv,pars=.BF_Fe2,twostage=TRUE,project=project)
}
BF_Ti <- function(A,T,M,project=FALSE){
    uv <- alr(cbind(T,A,M))
    BF_helper(uv,pars=.BF_Ti,project=project)
}
BF_helper <- function(uv,pars,twostage=TRUE,project=FALSE){
    if (twostage){
        xica <- pars$p[1]
        yica <- pars$p[2]
        xith <- pars$p[3]
        yith <- pars$p[4]
        xib <- (xica+xith)/2
        yib <- (yica+yith)/2
        b1 <- pars$p[5]
        b2 <- pars$p[6]
        ns <- nrow(uv)
        if (project){
            out <- matrix(NA,nrow=ns,ncol=2)
        } else {
            out <- rep(NA,ns)
        }
        for (i in 1:ns){
            xy <- uv[i,]
            xypb <- nearestpoint(xy,c(xib,yib,b1,b2))
            if (b1<0){
                if (xypb[1]<xib) {
                    b <- b1
                } else {
                    b <- b2
                }
            } else {
                if (xypb[2]>yib){
                    b <- b1
                } else {
                    b <- b2
                }
            }
            xyca <- point2line(xy,aa=b,bb=-1,cc=yica-b*xica)
            xyth <- point2line(xy,aa=b,bb=-1,cc=yith-b*xith)
            xyb <- point2line(xy,aa=b,bb=-1,cc=yib-b*xib)
            xyip <- point2line(c(xica,yica),aa=b,bb=-1,cc=yith-b*xith)
            D <- sqrt(sum((c(xica,yica)-xyip)^2))/2
            dca <- sqrt(sum(xyca-xy)^2)
            dth <- sqrt(sum(xyth-xy)^2)
            d <- sqrt(sum(xyb-xy)^2)
            if (project){
                out[i,] <- xypb
            } else if (dca<dth){
                out[i] <- d/D
            } else {
                out[i] <- -d/D
            }
        }
    } else {
        out <- pars$p[1]*atan((uv[,2]-pars$xy0[2])/(uv[,1]-pars$xy0[1]))+pars$p[2]
    }
    out
}
