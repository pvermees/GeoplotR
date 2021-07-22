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

BF <- function(A,F,M,T,twostage=TRUE){
    if (!missing(F)){
        out <- BF_Fe(A,F,M,twostage=twostage)
    } else if (!missing(T)){
        out <- BF_Ti(A,T,M)
    } else {
        stop('The Bowen-Fenner index requires either F or T.')
    }
    out
}
BF_Fe <- function(A,F,M,twostage=TRUE,project=FALSE){    
    uv <- alr(cbind(F,A,M))
    p <- c(-0.8,-0.3,-1.45,-1,-6,-0.6)
    BF_helper(uv,p=p,twostage=twostage,project=project)
}
BF_Ti <- function(A,T,M,project=FALSE){
    uv <- alr(cbind(T,A,M))
    p <- c(1,1.4,0,0.65,2.5,-0.35)
    BF_helper(uv,p=p,project=project)
}
BF_helper <- function(uv,p=p,twostage=FALSE,project=FALSE){
    if (twostage){
        xica <- p[1]
        yica <- p[2]
        xith <- p[3]
        yith <- p[4]
        xib <- (xica+xith)/2
        yib <- (yica+yith)/2
        b1 <- p[5]
        b2 <- p[6]
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
        out <- 9.57 * atan((uv[,2]-2.12)/(uv[,1]+3.84)) + 6.98
    }
    out
}

#' @title A-F-M
#' @description A-F-M diagram
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param F vector with (FeO + FeO\eqn{_2}O\eqn{_3}) concentrations
#'     (in wt\%)
#' @param M vector with MgO concentrations (in wt\%)
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#' @param twostage logical. If \code{TRUE}, applies the two-stage
#'     magma evolution model of Vermeesch and Pease (2021). Otherwise
#'     applies the single stage model.
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
#' @param xpd clip the plot region? See \code{par} for further details.
#' @param bw the smoothing bandwidth to be used for the kernel density
#'     estimate. See \code{density}.
#' @param dlty line type of the decision boundary.
#' @param dlwd line width of the decision boundary.
#' @param dcol colour of the decision boundary.
#' @param padding fractional measure of distance between the data and
#'     the kernel density estimate.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return the Bowen-Fenner indices of the samples, where positive
#'     values indicate calc-alkaline and negative numbers tholeiitic
#'     compositions (Vermeesch and Pease, 2021).
#' @references Irvine, T. N. and Baragar, W. A guide to the chemical
#'     classification of the common volcanic rocks. Canadian Journal
#'     of Earth Sciences, 8(5):523--548, 1971.
#'
#' Kuno, H. Differentiation of basalt magmas. Basalts: The Poldervaart
#'     treatise on rocks of basaltic composition, pages 623--688, 1968.
#'
#' Vermeesch, P. and Pease, V. A genetic classification of the
#' tholeiitic and calc-alkaline magma series, Geochemical Perspective
#' Letters (in review).
#' 
#' @examples
#' data(cath,package='GeoplotR')
#' oldpar <- par(mfrow=c(2,2))
#' on.exit(par(oldpar))
#'
#' A <- cath[,'Na2O']+cath[,'K2O']
#' F <- cath[,'FeOT']
#' M <- cath[,'MgO']
#' tern <- c(TRUE,FALSE,TRUE,FALSE)
#' stag <- c(TRUE,TRUE,FALSE,FALSE)
#' for (i in 1:4){
#'     AFM(A,F,M,ternary=tern[i],decision=TRUE,twostage=stag[i],
#'         pch=21,bw=0.2,padding=0.15,bg=cath[,'affinity'])
#' }
#'
#' dev.new()
#' CA <- (cath$affinity=='ca')
#' TH <- (cath$affinity=='th')
#' fitCA <- AFM(A[CA],F[CA],M[CA],plot=FALSE)
#' fitTH <- AFM(A[TH],F[TH],M[TH],plot=FALSE)
#' d1 <- density(fitCA)
#' d2 <- density(fitTH)
#' matplot(cbind(d1$x,d2$x),cbind(d1$y,d2$y),type='l',lty=1)
#' legend('topright',legend=c('Cascades','Iceland'),lty=1)
#' 
#' @export
AFM <- function(A,F,M,ternary=TRUE,twostage=TRUE,plot=TRUE,kde=TRUE,
                decision=TRUE,bty='n',asp=1,xpd=FALSE,bw="nrd0",
                dlty=2,dlwd=1.5,dcol='blue',padding=0.15,...){
    miss <- (missing(A)|missing(F)|missing(M))
    if (miss){
        out <- NULL
        kde <- FALSE
    } else {
        out <- BF(A,F,M,twostage=twostage)
    }
    if (plot){
        if (miss | ternary){
            minu <- -10
            maxu <- 10
        } else {
            uv <- alr(cbind(F,A,M))
            if (twostage){
                uvp <- BF_Fe(A,F,M,twostage=twostage,project=TRUE)
                minu <- min(uvp[,1])
                maxu <- max(uvp[,1])+padding*(max(uvp[,1])+1.125)
            } else {
                xy0 <- c(-3.84,2.12)
                r <- sqrt(rowSums(sweep(uv,2,xy0,'-')^2))
                R <- max(r)*(1+padding)
                minu <- min(uv[,1])
                maxu <- xy0[1]+R*cos(-6.98/9.57)
            }
            minv <- min(uv[,2])
            maxv <- max(uv[,2])
        }
        if (decision){
            nuv <- 100
            if (twostage){
                u1 <- seq(from=min(-1.125,minu),to=-1.125,length.out=nuv)
                u2 <- seq(from=-1.125,to=max(-1.125,maxu),length.out=nuv)
                v1 <- -0.65-6*(u1+1.125)
                v2 <- -0.65-0.6*(u2+1.125)
                uvd <- cbind(c(u1,u2),c(v1,v2))
            } else {
                u <- seq(from=-3.84,to=max(-3.84,maxu),length.out=nuv)
                v <- 2.12 + tan(-6.98/9.57)*(u+3.84)
                uvd <- cbind(u,v)
            }
        }
        if (ternary){
            ternaryplot(labels=c('F','A','M'),xpd=xpd)
            if (!miss){
                xyz <- cbind(F=F,A=A,M=M)
                graphics::points(xyz2xy(xyz),...)
            }
            if (decision){
                graphics::lines(xyz2xy(alr(uvd,inverse=TRUE)),
                                lty=dlty,lwd=dlwd,col=dcol)
            }
        } else {
            if (kde){
                dens <- stats::density(out,bw=bw)
                if (twostage){
                    alpha <- atan(0.6)
                    # 1. density
                    d <- 0.5*(1.87-0.78)/sqrt(0.6^2+1)
                    D <- d*dens$x
                    xed <- maxu + D*sin(alpha)
                    yed <- (-0.65-0.6*(maxu+1.125)) + D*cos(alpha)
                    r <- (maxu-minu)/cos(alpha)
                    normdens <- r*padding*dens$y/max(dens$y)
                    xd <- xed + normdens*cos(alpha)
                    yd <- yed - normdens*sin(alpha)
                    # 2. ticks
                    ticks <- pretty(dens$x)
                    Dt <- d*ticks
                    xt0 <- maxu + Dt*sin(alpha)
                    yt0 <- (-0.65-0.6*(maxu+1.125)) + Dt*cos(alpha)
                    xt1 <- xt0 - r*0.02*cos(alpha)
                    yt1 <- yt0 + r*0.02*sin(alpha)
                    # 3. axis
                    xr <- range(xt0)
                    yr <- range(yt0)
                } else {
                    # 1. density
                    Rd <- R*padding*dens$y/max(dens$y)
                    alpha <- (dens$x-6.98)/9.57
                    xd <- xy0[1]+(R+Rd)*cos(alpha)
                    yd <- xy0[2]+(R+Rd)*sin(alpha)
                    # 2. ticks
                    ticks <- pretty(dens$x)
                    beta <- (ticks-6.98)/9.57
                    xt0 <- xy0[1]+R*cos(beta)
                    yt0 <- xy0[2]+R*sin(beta)                    
                    xt1 <- xy0[1]+R*0.98*cos(beta)
                    yt1 <- xy0[2]+R*0.98*sin(beta)
                    # 3. arc
                    gamma <- seq(from=min(beta),to=max(beta),length.out=50)
                    xr <- xy0[1]+R*cos(gamma)
                    yr <- xy0[2]+R*sin(gamma)
                }
                maxu <- max(xd)
                minv <- min(yd)
            }
            if (!miss){ 
                graphics::plot(x=c(minu,maxu),y=c(minv,maxv),type='n',
                               xlab='ln(A/F)',ylab='ln(M/F)',
                               bty=bty,asp=asp,xpd=xpd,...)
                graphics::points(uv,...)
            } else {
                graphics::plot(c(-3,2),c(-3,1),type='n',xlab='ln(A/F)',
                               ylab='ln(M/F)',bty=bty,asp=asp,xpd=xpd,...)
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
