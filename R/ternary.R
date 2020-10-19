ternaryplot <- function(f=rep(1,3),labels=c('X','Y','Z')){
    corners <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
    xy <- xyz2xy(corners)
    graphics::plot(xy,type='l',asp=1,axes=FALSE,
                   ann=FALSE,bty='n')
    position <- c(3,1,1)
    for (i in 1:3){
        if (f[i]==1) lab <- labels[i]
        else lab <- paste0(f[i],'x',labels[i])
        graphics::text(xy[i,,drop=FALSE],labels=lab,pos=position[i])
    }
}
ternarylines <- function(uv,...){
    xyz <- alr(uv,inverse=TRUE)
    graphics::lines(xyz2xy(xyz),...)
}
ternarypoints <- function(uv,...){
    xyz <- alr(uv,inverse=TRUE)
    graphics::points(xyz2xy(xyz),...)
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
