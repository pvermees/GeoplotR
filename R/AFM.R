#' @title A-F-M
#' @description A-F-M diagram
#' @param A vector with (Na\eqn{_2}O+K\eqn{_2}O) concentrations (in
#'     wt\%)
#' @param F vector with (FeO + FeO\eqn{_2}O\eqn{_3}) concentrations
#'     (in wt\%)
#' @param M vector with MgO concentrations (in wt\%)
#' @param genetic logical. If \code{TRUE}, uses the logratio fit of
#'     Vermeesch and Pease (2021). Otherwise uses the empirical fit of
#'     Irvine and Baragar (1971).
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#' @param twostage logical. If \code{TRUE}, applies the two-stage
#'     magma evolution model of Vermeesch and Pease (2021). Otherwise
#'     applies the single stage model. Only used if
#'     \code{genetic=TRUE}.
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
#' @param pch plotting character. See the documentation of
#'     \code{points()} for details.
#' @param bg background colour for the plot symbols. If \code{NULL},
#'     uses the assigned classes.
#' @param dlwd line width of the decision boundary.
#' @param dcol colour of the decision boundary.
#' @param padding fractional measure of distance between the data and
#'     the kernel density estimate.
#' @param xlim the x axis limits of the plot (see \code{plot.default}
#'     for details).
#' @param ylim the y axis limits of the plot (see \code{plot.default}
#'     for details).
#' @param show.labels logical. If \code{TRUE}, labels the tholeiitic
#'     and calc-alkaline fields.
#' @param short logical. If \code{TRUE}, uses abbreviated labels (only
#'     used if \code{show.labels} is \code{TRUE}).
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
#' Irvine, T.N. and Baragar, W.R.A., 1971. A guide to the chemical
#' classification of the common volcanic rocks. Canadian journal of
#' earth sciences, 8(5), pp.523-548.
#'
#' Vermeesch, P. and Pease, V. A genetic classification of the
#' tholeiitic and calc-alkaline magma series, Geochemical Perspective
#' Letters.
#' 
#' @examples
#' data(cath,package='GeoplotR')
#'
#' A <- cath[,'Na2O']+cath[,'K2O']
#' F <- cath[,'FeOT']
#' M <- cath[,'MgO']
#' 
#' # 1. Irvine and Baragar (1971):
#' 
#' AFM(A=A,F=F,M=M,genetic=FALSE,
#'     bg=cath[,'affinity'],dcol='blue',pch=21)
#' 
#' # 2. Vermeesch and Pease (2021):
#' 
#' oldpar <- par(mfrow=c(2,2),mar=c(4,4,0,0))
#'
#' tern <- c(TRUE,FALSE,TRUE,FALSE)
#' stag <- c(TRUE,TRUE,FALSE,FALSE)
#' for (i in 1:4){
#'    AFM(A,F,M,ternary=tern[i],decision=TRUE,twostage=stag[i],
#'        pch=21,bw=0.2,padding=0.15,bg=cath[,'affinity'])
#' }
#'
#' par(oldpar)
#' 
#' @export
AFM <- function(A,F,M,genetic=TRUE,ternary=TRUE,twostage=TRUE,
                kde=TRUE,decision=TRUE,bty='n',asp=1,
                xpd=FALSE,bw="nrd0",pch=21,bg=NULL,dlwd=1.5,
                dcol='blue',padding=0.15,xlim=NULL,ylim=NULL,
                show.labels=FALSE,short=TRUE,...){
    miss <- (missing(A)|missing(F)|missing(M))
    if (miss){
        A <- NULL
        F <- NULL
        M <- NULL
        out <- NULL
        kde <- FALSE
    } else if (genetic){
        uv <- alr(cbind(F,A,M))
        out <- BF(A=A,F=F,M=M,twostage=twostage)
        if (is.null(bg)) bg <- ((out>0)+1)
    }
    if (genetic){
        if (miss | ternary){
            minu <- -10
            maxu <- 10
        } else {
            if (twostage){
                uvp <- BF_Fe2(A,F,M,project=TRUE)
                minu <- min(uvp[,1])
                maxu <- max(uvp[,1])+padding*(max(uvp[,1])-.BF_Fe2$xyi[1])
            } else {
                r <- sqrt(rowSums(sweep(uv,2,.BF_Fe1$xy0,'-')^2))
                R <- max(r)*(1+padding)
                minu <- min(uv[,1])
                maxu <- .BF_Fe1$xy0[1]+R*cos(-.BF_Fe1$p[2]/.BF_Fe1$p[1])
            }
            minv <- min(uv[,2])
            maxv <- max(uv[,2])
        }
        if (decision){
            nuv <- 100
            if (twostage){
                u1 <- seq(from=min(.BF_Fe2$xyi[1],minu),
                          to=.BF_Fe2$xyi[1],length.out=nuv)
                u2 <- seq(from=.BF_Fe2$xyi[1],
                          to=max(.BF_Fe2$xyi[1],maxu),length.out=nuv)
                v1 <- .BF_Fe2$xyi[2]+.BF_Fe2$b[1]*(u1-.BF_Fe2$xyi[1])
                v2 <- .BF_Fe2$xyi[2]+.BF_Fe2$b[2]*(u2-.BF_Fe2$xyi[1])
                uvd <- cbind(c(u1,u2),c(v1,v2))
            } else {
                u <- seq(from=.BF_Fe1$xy0[1],
                         to=max(.BF_Fe1$xy0[1],maxu),length.out=nuv)
                v <- .BF_Fe1$xy0[2] +
                    tan(-.BF_Fe1$p[2]/.BF_Fe1$p[1])*(u-.BF_Fe1$xy0[1])
                uvd <- cbind(u,v)
            }
        }
        if (ternary){
            ternaryplot(labels=c('F','A','M'),xpd=xpd)
            if (!miss){
                ternarypoints(uv,pch=pch,bg=bg,...)
            }
            if (decision){
                ternarylines(uvd,lty=1,lwd=dlwd,col=dcol)
            }
        } else {
            if (kde){
                dens <- stats::density(out,bw=bw)
                if (twostage){
                    alpha <- atan(-.BF_Fe2$b[2])
                    # 1. density
                    D <- .BF_Fe2$d*dens$x
                    xed <- maxu + D*sin(alpha)
                    yed <- D*cos(alpha) + .BF_Fe2$xyi[2] +
                        .BF_Fe2$b[2]*(maxu-.BF_Fe2$xyi[1])
                    r <- (maxu-minu)/cos(alpha)
                    normdens <- r*padding*dens$y/max(dens$y)
                    xd <- xed + normdens*cos(alpha)
                    yd <- yed - normdens*sin(alpha)
                    # 2. ticks
                    ticks <- pretty(dens$x)
                    Dt <- .BF_Fe2$d*ticks
                    xt0 <- maxu + Dt*sin(alpha)
                    yt0 <- Dt*cos(alpha) + .BF_Fe2$xyi[2] +
                        .BF_Fe2$b[2]*(maxu-.BF_Fe2$xyi[1])
                    xt1 <- xt0 - r*0.02*cos(alpha)
                    yt1 <- yt0 + r*0.02*sin(alpha)
                    # 3. axis
                    xr <- range(xt0)
                    yr <- range(yt0)
                } else {
                    # 1. density
                    Rd <- R*padding*dens$y/max(dens$y)
                    alpha <- (dens$x-.BF_Fe1$p[2])/.BF_Fe1$p[1]
                    xd <- .BF_Fe1$xy0[1]+(R+Rd)*cos(alpha)
                    yd <- .BF_Fe1$xy0[2]+(R+Rd)*sin(alpha)
                    # 2. ticks
                    ticks <- pretty(dens$x)
                    beta <- (ticks-.BF_Fe1$p[2])/.BF_Fe1$p[1]
                    xt0 <- .BF_Fe1$xy0[1]+R*cos(beta)
                    yt0 <- .BF_Fe1$xy0[2]+R*sin(beta)                    
                    xt1 <- .BF_Fe1$xy0[1]+R*0.98*cos(beta)
                    yt1 <- .BF_Fe1$xy0[2]+R*0.98*sin(beta)
                    # 3. arc
                    gamma <- seq(from=min(beta),to=max(beta),length.out=50)
                    xr <- .BF_Fe1$xy0[1]+R*cos(gamma)
                    yr <- .BF_Fe1$xy0[2]+R*sin(gamma)
                }
                maxu <- max(xt0)
                minv <- min(yt0)
            }
            if (!miss){
                graphics::plot(x=c(minu,maxu),y=c(minv,maxv),type='n',
                               xlab='ln(A/F)',ylab='ln(M/F)',
                               bty=bty,asp=asp,xpd=xpd,xlim=xlim,ylim=ylim)
                graphics::points(uv,pch=pch,bg=bg,...)
            } else {
                graphics::plot(c(-3,2),c(-3,1),type='n',xlab='ln(A/F)',
                               ylab='ln(M/F)',bty=bty,asp=asp,xpd=xpd,
                               xlim=xlim,ylim=ylim)
            }
            if (decision){
                graphics::lines(uvd,lty=2,lwd=dlwd,col=dcol)
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
        if (decision){
            plotlabels(diagram='AFM',ternary=ternary,
                       show.labels=show.labels,short=short)
        }
    } else {
        out <- xyzplot(json=.AFM,X=F,Y=A,Z=M,labels=c('F','A','M'),
                       dlwd=dlwd,dcol=dcol,pch=pch,bg=bg,
                       show.labels=show.labels,short=short,...)
    }
    invisible(out)
}
