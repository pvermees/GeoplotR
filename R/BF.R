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
        if (twostage){
            out <- BF_Fe2(A,F,M)
        } else {
            out <- BF_Fe1(A,F,M)
        }
    } else if (!missing(T)){
        out <- BF_Ti(A,T,M)
    } else {
        stop('The Bowen-Fenner index requires either F or T.')
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
BF_helper <- function(uv,pars,twostage=FALSE,project=FALSE){
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
