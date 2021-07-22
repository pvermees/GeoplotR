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
