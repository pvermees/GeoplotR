DA <- function(uv,da,D2=FALSE,ternary=FALSE,f=rep(1,3),tot=1,
               xlab=ifelse(D2,'1-X-Y','X'),
               ylab=ifelse(D2,'X','Y'),
               zlab=ifelse(D2,'Y','Z'),
               pch=21,bg=NULL,xlim=NULL,ylim=NULL,...){
    UV <- stats::na.omit(data.frame(u=uv[,1],v=uv[,2]))
    out <- DApredict(da$fit,UV)
    if (is.null(bg)) bg <- out$class
    if (ternary){
        p <- graphics::par(oma=rep(0,4),mar=rep(1,4),xpd=NA)
        if (D2) xyzlab <- c(paste0('1-',xlab,'-',ylab),xlab,ylab)
        else xyzlab <- c(xlab,ylab,zlab)
        ternaryplot(f=f,xyzlab=xyzlab)
        fcorr <- log(f[-1])-log(f[1])
        for (cont in da$contours){
            fcont <- sweep(cont,2,fcorr,'+')
            ternarylines(fcont)
        }
        fdat <- sweep(UV,2,fcorr,'+')
        ternarypoints(fdat,bg=bg,pch=pch,...)
        graphics::par(p)
    } else {
        if (D2){
            xlab <- xlab
            ylab <- ylab
            xy <- alr(UV,inverse=TRUE,tot=tot)[,-1]
        } else {
            den <- xlab
            xlab <- paste0('ln[',ylab,'/',den,']')
            ylab <- paste0('ln[',zlab,'/',den,']')
            xy <- UV
        }
        graphics::plot(xy,type='n',xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
        for (cont in da$contours){
            if (D2) graphics::lines(alr(cont,inverse=TRUE,tot=tot)[,-1])
            else graphics::lines(cont)
        }

        graphics::points(xy[,1:2],pch=pch,bg=bg,...)
    }
    invisible(out)
}

construct_DA <- function(X,Y,Z,quadratic=FALSE,plot=FALSE){
    out <- list()     
    x <- get_training_data(X)
    y <- get_training_data(Y)
    if (missing(Z)){
        out$ndim <- 2
        uv <- alr(cbind(1e6-x-y,x,y))
    } else {
        out$ndim <- 3
        uv <- alr(cbind(x,y,get_training_data(Z)))
    }
    u <- uv[,1]
    v <- uv[,2]
    dat <- data.frame(AFFINITY=training[,'AFFINITY'],u=u,v=v)
    nn <- 5000
    padding <- 4
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    vgrid <- seq(from=min(v,na.rm=TRUE)-padding,
                 to=max(v,na.rm=TRUE)+padding,length.out=nn)
    uvgrid <- expand.grid(ugrid,vgrid)
    if (quadratic){
        out$fit <- MASS::qda(AFFINITY ~ ., data=dat, na.action='na.omit')
        nt <- 750
    } else {
        out$fit <- MASS::lda(AFFINITY ~ ., data=dat, na.action='na.omit')
        nt <- 250
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
