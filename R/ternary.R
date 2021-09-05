ternaryplot <- function(labels=c('X','Y','Z'),f=rep(1,3),...){
    corners <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
    xy <- xyz2xy(corners)
    graphics::plot(xy,type='l',asp=1,axes=FALSE,
                   ann=FALSE,bty='n',...)
    position <- c(3,1,1)
    for (i in 1:3){
        if (f[i]==1) lab <- labels[i]
        else lab <- paste0(f[i],'x',labels[i])
        graphics::text(xy[i,,drop=FALSE],labels=lab,pos=position[i],...)
    }
}
ternarypoints <- function(uv,f=rep(1,3),...){
    ternaryhelper(uv=uv,type='p',f=f,...)
}
ternarylines <- function(uv,f=rep(1,3),...){
    ternaryhelper(uv=uv,type='l',f=f,...)
}
ternarytext <- function(uv,f=rep(1,3),labels=seq_along(uv[,1]),...){
    ternaryhelper(uv=uv,type='t',f=f,labels=labels,...)
}
ternaryhelper <- function(uv,type='p',f=rep(1,3),
                          labels=seq_along(uv[,1]),...){
    uvt <- uv - log(f[1])
    uvt <- sweep(uvt,2,log(f[2:3]),'+')
    xyz <- alr(uvt,inverse=TRUE)
    xy <- xyz2xy(xyz)
    if (type=='p'){
        graphics::points(xy,...)
    } else if (type=='l'){
        graphics::lines(xy,...)
    } else if (type=='t'){
        graphics::text(x=xy[,1],y=xy[,2],labels=labels,...)
    } else {
        stop('Invalid value for argument type in ternaryhelper().')
    }
}

xyz2xy <- function(xyz){
    if (class(xyz)%in%c('matrix','data.frame')){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    return(xy)
}
