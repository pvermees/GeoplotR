DA <- function(uv,da,ternary=FALSE,f=rep(1,3),
               xyzlab=c('X','Y','Z'),pch=21,bg=NULL,...){
    dat <- stats::na.omit(data.frame(u=uv[,1],v=uv[,2]))
    out <- DApredict(da$fit,dat)
    if (is.null(bg)) bg <- out$class
    if (ternary){
        p <- graphics::par(oma=rep(0,4),mar=rep(1,4),xpd=NA)
        ternaryplot(f=f,xyzlab=xyzlab)
        fcorr <- log(f[-1])-log(f[1])
        for (cont in da$contours){
            fcont <- sweep(cont,2,fcorr,'+')
            ternarylines(fcont)
        }
        fdat <- sweep(dat,2,fcorr,'+')
        ternarypoints(fdat,bg=bg,pch=pch,...)
        graphics::par(p)
    } else {
        xy <- do.call("rbind",da$contours)
        xlab <- paste0('ln[',xyzlab[2],'/',xyzlab[1],']')
        ylab <- paste0('ln[',xyzlab[3],'/',xyzlab[1],']')
        graphics::plot(xy,type='n',xlab=xlab,ylab=ylab)
        for (cont in da$contours){
            graphics::lines(cont)
        }
        graphics::points(x=dat[,1],y=dat[,2],pch=pch,bg=bg,...)
    }
    invisible(out)
}

DA2D <- function(uv,da,ternary=FALSE,f=rep(1,3),
                 xyzlab=c('1-X-Y','X','Y'),pch=21,bg=NULL,...){
    complete <- complete.cases(uv)
    dat <- stats::na.omit(data.frame(u=uv[,1],v=uv[,2]))
    out <- DApredict(da$fit,dat)
    if (is.null(bg)) bg <- out$class
    if (ternary){
        p <- graphics::par(oma=rep(0,4),mar=rep(1,4),xpd=NA)
        ternaryplot(f=f,xyzlab=xyzlab)
        fcorr <- log(f[-1])-log(f[1])
        for (cont in da$contours){
            fcont <- sweep(cont,2,fcorr,'+')
            ternarylines(fcont)
        }
        fdat <- sweep(dat,2,fcorr,'+')
        ternarypoints(fdat[complete,],bg=bg,pch=pch,...)
        graphics::par(p)
    } else {
        dlines <- do.call("rbind",da$contours)
        xlab <- xyzlab[2]
        ylab <- xyzlab[3]
        xy <- alr(uv,inverse=TRUE)[complete,-1]
        graphics::plot(xy,type='n',xlab=xlab,ylab=ylab)
        for (cont in da$contours){
            graphics::lines(alr(cont,inverse=TRUE)[,-1])
        }
        graphics::points(x=xy[,1],y=xy[,2],pch=pch,bg=bg,...)
    }
    invisible(out)
}

# units is only used if Z is missing
construct_DA <- function(X,Y,Z,quadratic=FALSE,plot=FALSE,units='ppm'){
    out <- list()
    AFFINITY <- training[,'AFFINITY']
    x <- get_training_data(X)
    y <- get_training_data(Y)
    if (missing(Z)){
        out$ndim <- 2
        if (identical(units,'ppm')){
            S <- 1e6
        } else {
            S <- 100
        }
        uv <- alr(cbind(S-x-y,x,y))
    } else {
        out$ndim <- 3
        uv <- alr(cbind(x,y,get_training_data(Z)))
    }
    nn <- 3000
    u <- uv[,1]
    v <- uv[,2]
    padding <- 4
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    vgrid <- seq(from=min(v,na.rm=TRUE)-padding,
                 to=max(v,na.rm=TRUE)+padding,length.out=nn)
    uvgrid <- expand.grid(ugrid,vgrid)
    if (quadratic){
        out$fit <- MASS::qda(AFFINITY ~ u + v,na.action='na.omit')
        nt <- 500
    } else {
        out$fit <- MASS::lda(AFFINITY ~ u + v,na.action='na.omit')
        nt <- 50
    }
    out$contours <- list()
    pr <- DApredict(fit=out$fit,dat=data.frame(u=uvgrid[,1],v=uvgrid[,2]))
    z <- matrix(as.numeric(pr$class),nrow=nn,ncol=nn)
    if (plot) plot(ugrid,vgrid,type='n')
    contours <- grDevices::contourLines(ugrid,vgrid,z,levels=c(1.5,2.5))
    for (i in seq_along(contours)){
        nx <- length(contours[[i]]$x)
        j <- seq(from=1,to=nx,length.out=nt)
        x <- contours[[i]]$x[j]
        y <- contours[[i]]$y[j]
        out$contours[[i]] <- cbind(x,y)
        if (plot) graphics::lines(x,y)
    }
    invisible(out)
}

LDApredict <- utils::getFromNamespace("predict.lda", "MASS")
QDApredict <- utils::getFromNamespace("predict.qda", "MASS")
DApredict <- function(fit,dat){
    if (class(fit)%in%'lda'){
        out <- LDApredict(fit,newdata=dat)
    } else {
        out <- QDApredict(fit,newdata=dat)
    }
}
