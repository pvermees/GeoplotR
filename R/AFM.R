# aa*x + bb*y + cc = 0
point2line <- function(xy,aa,bb,cc){
    x <- xy[1]
    y <- xy[2]
    x0 <- (bb*(bb*x-aa*y)-aa*cc)/(aa^2+bb^2)
    y0 <- (aa*(aa*y-bb*x)-bb*cc)/(aa^2+bb^2)
    c(x0,y0)
}
lAMFi <- function(p){
    out <- list()
    out$a1 <- p[1]-p[2]*p[3]
    out$b1 <- p[3]
    out$a2 <- p[1]+(1-p[4])*p[5]-p[2]*p[4]
    out$b2 <- p[4]
    out$xy0 <- p[1:2]
    out$xyi <- c((out$a2-out$a1)/(out$b1-out$b2),
    (out$b1*out$a2-out$a1*out$b2)/(out$b1-out$b2))
    xyp <- point2line(out$xyi,aa=out$b2,bb=-1,cc=out$xy0[2]-out$b2*out$xy0[1])
    out$xy1 <- (out$xyi+xyp)/2
    out
}
# xy = 2-column matrix
nearestpoints <- function(xy,par,th=TRUE){
    # xy = 2-element vector
    nearestpoint <- function(xy,par,th=TRUE){
        p <- lAMFi(par)
        if (th){
            xyp1 <- point2line(xy,p$b1,-1,p$a1)
            xyp2 <- point2line(xy,p$b2,-1,p$a2)
            if (xyp1[2]>p$xy0[2]){ # above initial ratios
                out <- p$xy0
            } else {
                offside1 <- (xyp1[2]<p$xyi[2])
                offside2 <- (xyp2[1]<p$xyi[1])
                if (offside1){
                    if (offside2) out <- p$xyi
                    else out <- xyp2
                } else if (offside2){
                    out <- xyp1
                } else {
                    d1 <- sqrt(sum((xy-xyp1)^2))
                    d2 <- sqrt(sum((xy-xyp2)^2))
                    if (d1<d2) out <- xyp1
                    else out <- xyp2
                }
            }
        } else {
            out <- point2line(xy,par[4],-1,par[1]-par[4]*par[2])
            if (out[1]<p$xy0[1]) out <- p$xy0
        }
        out
    }
    t(apply(xy,1,nearestpoint,par=par,th=th))
}
fitcath <- function(init,uvca,uvth){
    misfit <- function(p,xy,th=TRUE){
        xy0 <- nearestpoints(xy,p,th=th)
        d <- sqrt(rowSums((xy-xy0)^2))
        sum(d)
    }
    cathmisfit <- function(p,uvca,uvth){
        misfit(p,uvca,th=FALSE) + misfit(p,uvth,th=TRUE)
    }
    optim(init,cathmisfit,uvca=uvca,uvth=uvth)$par
}

construct_twostage <- function(cath,p=c(0.5,-2,-3,-0.5,-1)){
    FAM <- cbind(cath[,'FeOT'],cath[,'Na2O']+cath[,'K2O'],cath[,'MgO'])
    uv <- alr(FAM)    
    ith <- which(cath[,'affinity']=='th')
    ica <- which(cath[,'affinity']=='ca')
    fit <- fitcath(p,uvca=uv[ica,],uvth=uv[ith,])
    lAMFi(fit)
}
construct_BF <- function(cath){
    orthofit <- function(ab,uv){
        uvp <- t(apply(uv,1,point2line,aa=ab[2],bb=-1,cc=ab[1]))
        sqrt(sum((uv-uvp)^2))
    }
    FAM <- cbind(cath[,'FeOT'],cath[,'Na2O']+cath[,'K2O'],cath[,'MgO'])
    uv <- alr(FAM)
    ith <- which(cath[,'affinity']=='th')
    ica <- which(cath[,'affinity']=='ca')
    cafit <- optim(c(0,0),orthofit,uv=uv[ith,])$par
    thfit <- optim(c(0,0),orthofit,uv=uv[ica,])$par
    m0 <- tan((atan(cafit[2])+atan(thfit[2]))/2)
    dm <- tan((atan(cafit[2])-atan(thfit[2]))/2)
    x0 <- (cafit[1]-thfit[1])/(thfit[2]-cafit[2])
    y0 <- thfit[2]*x0 + thfit[1]
    list(th0=thfit[1],ca0=cafit[1],mth=thfit[2],
         mca=cafit[2],m0=m0,dm=dm,x0=x0,y0=y0)
}

#' @title Bowen-Fenner index
#' @description Bowen-Fenner index for tholeiitic and calc-alkaline
#'     magma series
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param F vector with (FeO + FeO\eqn{_2}O\eqn{_3}) concentrations
#'     (in wt\%)
#' @param M vector with MgO concentrations (in wt\%)
#' @return the Bowen-Fenner indices of the samples, where positive
#'     values indicate calc-alkaline and negative numbers tholeiitic
#'     compositions (Vermeesch and Pease, 2021).
#' @references Vermeesch, P. and Pease, V. A genetic classification of
#'     the tholeiitic and calc-alkaline magma series, Geochemical
#'     Perspective Letters (in review).
#' @examples
#' data(cath,package='GeoplotR')
#' bfi <- BF(A=cath$Na2O+cath$K2O,F=cath$FeOT,M=cath$MgO)
#' oldpar <- par(mfrow=c(2,1))
#' on.exit(par(oldpar))
#' hist(bfi[cath$affinity=='ca'])
#' hist(bfi[cath$affinity=='th'])
#' @export
BF <- function(A,F,M){
    m <- (log(M/F)-.BF$y0)/(log(A/F)-.BF$x0)
    (atan(m)-atan(.BF$m0))/atan(.BF$dm)
}

#' @title A-F-M
#' @description A-F-M diagram
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param F vector with (FeO + FeO\eqn{_2}O\eqn{_3}) concentrations
#'     (in wt\%)
#' @param M vector with MgO concentrations (in wt\%)
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#' @param radial logical. If \code{TRUE}, adds a kernel density
#'     estimate of the Bowen-Fenner indices on a radial scale. Not
#'     used if \code{ternary=FALSE}.
#' @param bty A character string that determines the type of `box'
#'     which is drawn about plots. See \code{par}.
#' @param asp the y/x aspect ratio, see \code{plot.window}.
#' @param bw the smoothing bandwidth to be used for the kernel density
#'     estimate. See \code{density}.
#' @param decision logical. If \code{TRUE}, adds Vermeesch and Pease's
#'     two-stage decision boundary between the calc-alkaline and
#'     tholeiitic magma series to the plot.
#' @param dlty line type of the decision boundary.
#' @param dlwd line width of the decision boundary.
#' @param dcol colour of the decision boundary.
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
#' Cascades <- (cath$affinity=='ca')
#' AFM(A=(cath$Na2O+cath$K2O)[Cascades],
#'     F=cath$FeOT[Cascades],M=cath$MgO[Cascades])
#' @export
AFM <- function(A,F,M,ternary=TRUE,radial=FALSE,bty='n',asp=1,xpd=FALSE,
                bw="nrd0",decision=TRUE,dlty=2,dlwd=1.5,dcol='blue',...){
    miss <- (missing(A)|missing(F)|missing(M))
    if (decision){
        nuv <- 20
        u1 <- seq(from=.twostage$xy0[1],to=.twostage$xy1[1],length.out=nuv)
        v1 <- seq(from=.twostage$xy0[2],to=.twostage$xy1[2],length.out=nuv)
        u2 <- seq(from=.twostage$xy1[1],to=10*diff(.twostage$xy1),length.out=nuv)
        v2 <- .twostage$xy1[2] + (u2-.twostage$xy1[1])*.twostage$b2
        uv <- cbind(c(u1,u2),c(v1,v2))
    }
    if (ternary){
        ternaryplot(labels=c('F','A','M'),xpd=xpd)
        if (!miss){
            xyz <- cbind(F=F,A=A,M=M)
            graphics::points(xyz2xy(xyz),...)
        }
        if (decision){
            graphics::lines(xyz2xy(alr(uv,inverse=TRUE)),
                            lty=dlty,lwd=dlwd,col=dcol)
        }
    } else {
        if (radial){
            radialBF(A,F,M,bty=bty,bw=bw,asp=asp,xpd=xpd,...)
        } else {
            if (!miss){
                graphics::plot(x=log(A/F),y=log(M/F),xlab='ln(A/F)',
                               ylab='ln(M/F)',bty=bty,asp=asp,xpd=xpd,...)
            } else {
                graphics::plot(c(-3,2),c(-3,1),type='n',xlab='ln(A/F)',
                               ylab='ln(M/F)',bty=bty,asp=asp,xpd=xpd,...)
            }
            if (decision){
                graphics::lines(uv,lty=dlty,lwd=dlwd,col=dcol)
            }
        }
    }
    if (miss){
        out <- NULL
    } else if (radial){
        out <- BF(A,F,M)
    } else {
        x <- log(A/F)
        y <- log(M/F)
        left <- (x < .twostage$xy1[1])
        right <- !left
        below <- function(x,y,a,b){ y < a + b*x }
        isTH <- ((left & below(x,y,.twostage$a1,.twostage$b1)) |
                 (right & below(x,y,.twostage$a2,.twostage$b2)))
        out <- rep('CA',length(A))
        out[isTH] <- "TH"
    }
    invisible(out)
}

radialBF <- function(A,F,M,bw="nrd0",asp=1,bty='n',...){
    x <- log(A/F)
    y <- log(M/F)
    m <- atan((y-.BF$y0)/(x-.BF$x0))
    dens <- stats::density(m,bw=bw)
    r <- 1.1*max(sqrt((x-.BF$x0)^2+(y-.BF$y0)^2))
    f <- .3
    X <- .BF$x0+(r+f*dens$y)*cos(dens$x)
    Y <- .BF$y0+(r+f*dens$y)*sin(dens$x)
    graphics::plot(x=c(.BF$x0,x,X),y=c(.BF$y0,y,Y),type='n',
                   xlab='ln(A/F)',ylab='ln(M/F)',
                   bty=bty,asp=asp)
    graphics::lines(X,Y,type='l')
    xrca <- .BF$x0+c(0,r*cos(atan(.BF$mca)))
    yrca <- .BF$y0+c(0,r*sin(atan(.BF$mca)))
    xr0 <- .BF$x0+c(0,r*cos(atan(.BF$m0)))
    yr0 <- .BF$y0+c(0,r*sin(atan(.BF$m0)))
    xrth <- .BF$x0+c(0,r*cos(atan(.BF$mth)))
    yrth <- .BF$y0+c(0,r*sin(atan(.BF$mth)))
    xr2ca <- .BF$x0+r*cos(2*atan(.BF$mca)-atan(.BF$m0))*c(1,0.98)
    yr2ca <- .BF$y0+r*sin(2*atan(.BF$mca)-atan(.BF$m0))*c(1,0.98)
    xr2th <- .BF$x0+r*cos(2*atan(.BF$mth)-atan(.BF$m0))*c(1,0.98)
    yr2th <- .BF$y0+r*sin(2*atan(.BF$mth)-atan(.BF$m0))*c(1,0.98)
    graphics::lines(xrca,yrca)
    graphics::lines(xr0,yr0,lty=2)
    graphics::lines(xrth,yrth)
    graphics::lines(xr2ca,yr2ca)
    graphics::lines(xr2th,yr2th)
    M <- seq(from=atan((Y[1]-.BF$y0)/(X[1]-.BF$x0)),
             to=atan((utils::tail(Y,1)-.BF$y0)/(utils::tail(X,1)-.BF$x0)),
             length.out=100)
    graphics::lines(x=.BF$x0+r*cos(M),y=.BF$y0+r*sin(M))
    graphics::text(xr2th[2],yr2th[2],labels='-2',pos=2)
    graphics::text(xr2ca[2],yr2ca[2],labels='+2',pos=2)
    graphics::points(x,y,...)
}
