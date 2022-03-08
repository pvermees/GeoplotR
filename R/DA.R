DA <- function(uv,da,ternary=TRUE,f=rep(1,3),
               xyzlab=c('X','Y','Z'),pch=21,bg=NULL,...){
    dat <- stats::na.omit(data.frame(u=uv[,1],v=uv[,2]))
    out <- DApredict(da$fit,dat)
    if (is.null(bg)) bg <- out$class
    else bg <- rep(1,nrow(uv))
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

construct_DA <- function(X,Y,Z,quadratic=FALSE,plot=FALSE){
    AFFINITY <- training[,'AFFINITY']
    x <- get_training_data(X)
    y <- get_training_data(Y)
    z <- get_training_data(Z)
    nn <- 3000
    uv <- alr(cbind(x,y,z))
    u <- uv[,1]
    v <- uv[,2]
    padding <- 4
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    vgrid <- seq(from=min(v,na.rm=TRUE)-padding,
                 to=max(v,na.rm=TRUE)+padding,length.out=nn)
    uvgrid <- expand.grid(ugrid,vgrid)
    out <- list()
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

construct_DA_2D <- function(X,Y,quadratic=FALSE,plot=FALSE){
    AFFINITY <- training[,'AFFINITY']
    x <- get_training_data(X)
    y <- get_training_data(Y)
    nn <- 1000
    u <- alr(cbind(x,y))
    padding <- 4
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    out <- list()
    if (quadratic){
        out$fit <- MASS::qda(AFFINITY ~ u,na.action='na.omit')
        nt <- 500
    } else {
        out$fit <- MASS::lda(AFFINITY ~ u,na.action='na.omit')
        nt <- 50
    }
    pr <- DApredict(fit=out$fit,dat=data.frame(u=ugrid))
    naff <- length(levels(pr$class))
    out$contours <- rep(NA,naff-1)
    for (i in 1:(naff-1)){
        aff <- pr$class[i]
        i <- tail(which(pr$class %in% aff),n=1)
        out$contours[i] <- mean(pr$x[i+(0:1)])
    }
    if (plot) matplot(ugrid,pr$posterior,type='l',lty=1,col='black')
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
