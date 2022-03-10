#' @title A-T-M
#' @description A-T-M diagram
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param T vector with TiO\eqn{_2} concentrations
#' @param M vector with MgO concentrations (in wt\%)
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#' @param plot logical. If \code{FALSE}, omits the graphical output
#'     and only returns the Bowen-Fenner indices of the sample.
#' @param kde logical. If \code{TRUE}, adds a kernel density estimate
#'     of the Bowen-Fenner indices on a radial scale. Not used if
#'     \code{ternary=FALSE}.
#' @param decision logical. If \code{TRUE}, adds the decision boundary
#'     between the tholeiitic and calc-alkaline magma series to the
#'     plot.
#' @param bty A character string that determines the type of `box'
#'     which is drawn about plots. See \code{par}.
#' @param asp the y/x aspect ratio, see \code{plot.window}.
#' @param xpd clip the plot region? See \code{par} for further
#'     details.
#' @param bw the smoothing bandwidth to be used for the kernel density
#'     estimate. See \code{density}.
#' @param dlty line type of the decision boundary.
#' @param dlwd line width of the decision boundary.
#' @param dcol colour of the decision boundary.
#' @param padding fractional measure of distance between the data and
#'     the kernel density estimate.
#' @param xlim the x axis limits of the plot (see \code{plot.default}
#'     for details).
#' @param ylim the y axis limits of the plot (see \code{plot.default}
#'     for details).
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return the Bowen-Fenner indices of the samples, where positive
#'     values indicate calc-alkaline and negative numbers tholeiitic
#'     compositions (Vermeesch and Pease, 2021).
#' @references Vermeesch, P. and Pease, V. A genetic classification of the
#' tholeiitic and calc-alkaline magma series, Geochemical Perspective
#' Letters.
#' 
#' @examples
#' data(cath,package='GeoplotR')
#' 
#' oldpar <- par(mfrow=c(2,2),mar=c(4,4,0,0))
#' on.exit(par(oldpar))
#'
#' A <- cath[,'Na2O']+cath[,'K2O']
#' T <- cath[,'TiO2']
#' M <- cath[,'MgO']
#' 
#' ATM(A,T,M,decision=FALSE,pch=21,bg='white',cex=0.5)
#'
#' ATM(A,T,M,ternary=TRUE,pch=21,bw=0.2,
#'     padding=0.15,bg=cath[,'affinity'])
#'
#' ATM(A,T,M,ternary=FALSE,decision=TRUE,kde=TRUE,
#'     pch=21,bw=0.2,padding=0.15,bg=cath[,'affinity'])
#'    
#' CA <- (cath$affinity=='ca')
#' TH <- (cath$affinity=='th')
#' fitCA <- ATM(A[CA],T[CA],M[CA],plot=FALSE)
#' fitTH <- ATM(A[TH],T[TH],M[TH],plot=FALSE)
#' d1 <- density(fitCA)
#' d2 <- density(fitTH)
#' matplot(cbind(d1$x,d2$x),cbind(d1$y,d2$y),
#'         type='l',lty=1,xlab='BF',ylab='')
#' legend('topleft',legend=c('Cascades','Iceland'),lty=1)
#'
#' @export
ATM <- function(A,T,M,ternary=TRUE,plot=TRUE,kde=TRUE,decision=TRUE,
                bty='n',asp=1,xpd=FALSE,bw="nrd0",pch=21,bg=NULL,
                dlty=2,dlwd=1.5,dcol='blue',padding=0.15,xlim=NULL,
                ylim=NULL,...){
    miss <- (missing(A)|missing(T)|missing(M))
    if (miss){
        out <- NULL
        kde <- FALSE
    } else {
        uv <- alr(cbind(T,A,M))
        out <- BF(A=A,T=T,M=M)
        if (is.null(bg)) bg <- ((out>0)+1)
    }
    if (plot){
        if (miss | ternary){
            minu <- -10
            maxu <- 10
        } else {
            uvp <- BF_Ti(A,T,M,project=TRUE)
            minu <- min(uvp[,1])
            maxu <- max(uvp[,1])+padding*(max(uvp[,1])-.BF_Ti$xyi[1])
            minv <- min(uv[,2])
            maxv <- max(uv[,2])
        }
        if (decision){
            nuv <- 100
            u1 <- seq(from=max(.BF_Ti$xyi[1],maxu),
                      to=.BF_Ti$xyi[1],length.out=nuv)
            u2 <- seq(from=.BF_Ti$xyi[1],
                      to=max(.BF_Ti$xyi[1],maxu),length.out=nuv)
            v1 <- .BF_Ti$xyi[2]+.BF_Ti$b[1]*(u1-.BF_Ti$xyi[1])
            v2 <- .BF_Ti$xyi[2]+.BF_Ti$b[2]*(u2-.BF_Ti$xyi[1])
            uvd <- cbind(c(u1,u2),c(v1,v2))
        }
        if (ternary){
            f <- c(5,1,1)
            ternaryplot(xyzlab=c('T','A','M'),f=f,xpd=xpd)
            if (!miss){
                ternarypoints(uv,f=f,pch=pch,bg=bg,...)
            }
            if (decision){
                ternarylines(uvd,f=f,lty=dlty,lwd=dlwd,col=dcol)
            }
        } else {
            if (kde){
                dens <- stats::density(out,bw=bw)
                alpha <- atan(-.BF_Ti$b[2])
                # 1. density
                D <- .BF_Ti$d*dens$x
                xed <- maxu + D*sin(alpha)
                yed <- D*cos(alpha) + .BF_Ti$xyi[2] +
                    .BF_Ti$b[2]*(maxu-.BF_Ti$xyi[1])
                r <- (maxu-minu)/cos(alpha)
                normdens <- r*padding*dens$y/max(dens$y)
                xd <- xed + normdens*cos(alpha)
                yd <- yed - normdens*sin(alpha)
                # 2. ticks
                ticks <- pretty(dens$x)
                Dt <- .BF_Ti$d*ticks
                xt0 <- maxu + Dt*sin(alpha)
                yt0 <- Dt*cos(alpha) + .BF_Ti$xyi[2] +
                    .BF_Ti$b[2]*(maxu-.BF_Ti$xyi[1])
                xt1 <- xt0 - r*0.02*cos(alpha)
                yt1 <- yt0 + r*0.02*sin(alpha)
                # 3. axis
                xr <- range(xt0)
                yr <- range(yt0)
                maxu <- xr[2]
                minv <- yr[1]
            }
            if (!miss){
                graphics::plot(x=c(minu,maxu),y=c(minv,maxv),type='n',
                               xlab='ln(A/T)',ylab='ln(M/T)',
                               bty=bty,asp=asp,xpd=xpd,xlim=xlim,ylim=ylim)
                graphics::points(uv,pch=pch,bg=bg,...)
            } else {
                graphics::plot(c(-1,3),c(-1,3),type='n',xlab='ln(A/T)',
                               ylab='ln(M/T)',bty=bty,asp=asp,xpd=xpd,
                               xlim=xlim,ylim=ylim)
            }
            if (decision){
                graphics::lines(uvd,lty=dlty,lwd=dlwd,col=dcol)
            }
            if (kde){
                graphics::lines(xr,yr)
                graphics::lines(xd,yd)
                graphics::matlines(x=rbind(xt0,xt1),
                                   y=rbind(yt0,yt1),
                                   lty=1,col='black')
                graphics::text(x=xt1,y=yt1,labels=ticks,pos=2)
            }
        }
    }
    invisible(out)
}
