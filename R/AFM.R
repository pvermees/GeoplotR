construct_BF <- function(cath){
    uv <- alr(cath[,-1])
    ith <- which(cath[,'affinity']=='th')
    ica <- which(cath[,'affinity']=='ca')
    th <- uv[ith,]
    ca <- uv[ica,]
    thfit <- stats::lm(th[,'M'] ~ th[,'A'])
    cafit <- stats::lm(ca[,'M'] ~ ca[,'A'])
    th0 <- thfit$coefficients[1]
    ca0 <- cafit$coefficients[1]
    mth <- thfit$coefficients[2]
    mca <- cafit$coefficients[2]
    m0 <- tan((atan(mca)+atan(mth))/2)
    dm <- tan((atan(mca)-atan(mth))/2)
    x0 <- (ca0-th0)/(mth-mca)
    y0 <- mth*x0 + th0
    list(th0=th0,ca0=ca0,mth=mth,mca=mca,m0=m0,dm=dm,x0=x0,y0=y0)
}

#' @title BF
#' @description Bowen-Fenner index
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
#' bfi <- BF(cath$A,cath$F,cath$M)
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
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return the Bowen-Fenner indices of the samples, where positive
#'     values indicate calc-alkaline and negative numbers tholeiitic
#'     compositions (Vermeesch and Pease, 2021).
#' @references Irvine, T. N. and Baragar, W. A guide to the chemical
#'     classification of the common volcanic rocks. Canadian Journal
#'     of Earth Sciences, 8(5):523–548, 1971.
#'
#' Kuno, H. Differentiation of basalt magmas. Basalts: The Poldervaart
#'     treatise on rocks of basaltic composition, pages 623–688, 1968.
#'
#' Vermeesch, P. and Pease, V. A genetic classification of the
#' tholeiitic and calc-alkaline magma series, Geochemical Perspective
#' Letters (in review).
#' 
#' @examples
#' data(cath,package='GeoplotR')
#' Cascades <- (cath$affinity=='ca')
#' AFM(cath$A[Cascades],cath$F[Cascades],cath$M[Cascades],
#'     ternary=FALSE,radial=TRUE)
#' @export
AFM <- function(A,F,M,ternary=TRUE,radial=FALSE,bty='n',asp=1,bw="nrd0",...){
    if (ternary){
        ternaryplot(labels=c('F','A','M'))
        xyz <- cbind(F=F,A=A,M=M)
        graphics::points(xyz2xy(xyz),...)
    } else {
        if (radial){
            radialBF(A,F,M,bty=bty,bw=bw,asp=asp,...)
        } else {
            graphics::plot(x=log(A/F),y=log(M/F),
                           xlab='ln(A/F)',ylab='ln(M/F)',
                           bty=bty,asp=asp,...)
        }
    }
    invisible(BF(A,F,M))
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
