AlFeTiMg <- function(Al,Fe,Ti,Mg,...){
    xyplot(x=Al+Fe,y=Ti+Mg,coords=.AlFeTiMg$coords,...)
}

AnAbOr <- function(An,Ab,Or,...){
    ternaryplot(labels=c('An','Ab','Or'))
    lnames <- names(.AnAbOr$polygons)
    coords <- .AnAbOr$polygons
    for (lname in lnames){
#        if (.AnAbOr$type[[lname]]=="dashed") lty <- 2
#        else if (.AnAbOr$type[[lname]]=="dotted") lty <- 3
#        else lty <- 1
        xyz <- matrix(unlist(coords[[lname]]),ncol=3,byrow=TRUE)
        graphics::lines(xyz2xy(xyz),lty=1)
    }
    lnames <- names(.AnAbOr$labels)
    coords <- .AnAbOr$labels
    angles <- .AnAbOr$angle
    for (lname in lnames){
        graphics::text(xyz2xy(coords[[lname]]),
                       labels=lname,srt=angles[[lname]])
    }
}

TiZrY_nominal <- function(Ti=NULL,Zr=NULL,Y=NULL,show.classes=TRUE,...){
    p <- graphics::par(oma=rep(0,4),mar=rep(1,4),xpd=NA)
    f <- c(0.01,1,3)
    ternaryplot(f=f,labels=c('Ti','Zr','Y'))
    lcoords <- .TiZrY_nominal$line_coords
    ltype <- .TiZrY_nominal$line_type
    for (lname in names(lcoords)){
        xyz <- matrix(unlist(lcoords[[lname]]),ncol=3,byrow=TRUE)
        graphics::lines(xyz2xy(xyz),lty=lty(ltype[[lname]]))
    }
    pcoords <- .TiZrY_nominal$polygons
    classes <- .TiZrY_nominal$labels
    pnames <- names(pcoords)
    if (any(is.null(Ti) | is.null(Zr) | is.null(Y))){
        out <- NULL
    } else {
        Z <- 3*Y
        Y <- Zr
        X <- Ti/100
        XYZ <- cbind(X,Y,Z)
        uv <- alr(XYZ,inverse=FALSE)
        XY <- xyz2xy(alr(uv,inverse=TRUE))
        ns <- nrow(XY)
        out <- rep(NA,ns)
        col <- rep(1,ns)
        for (i in 1:length(pcoords)){
            pname <- pnames[i]
            xyz <- matrix(unlist(pcoords[[pname]]),ncol=3,byrow=TRUE)
            matched <- apply(XY,MARGIN=1,FUN=in_polygon,pol=xyz2xy(xyz))
            out[matched] <- classes[[pname]]
            col[matched] <- i+1
        }
        ternarypoints(uv,col=col,...)
    }
    if (show.classes){
        for (i in 1:length(pcoords)){
            pname <- pnames[i]
            xyz <- matrix(unlist(.TiZrY_nominal$label_coords[[pname]]),ncol=3,byrow=TRUE)
            graphics::text(xyz2xy(xyz),labels=classes[[pname]],
                           srt=.TiZrY_nominal$angle[[pname]],pos=1)#,col=i+1)
        }
    }
    graphics::par(p)
    out
}
