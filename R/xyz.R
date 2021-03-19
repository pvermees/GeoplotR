AnAbOr <- function(An=NULL,Ab=NULL,Or=NULL,show.labels=TRUE,...){
    invisible(xyzplot(json=.AnAbOr,X=An,Y=Ab,Z=Or,labels=c('An','Ab','Or'),
                      show.labels=show.labels,...))
}

TiZrY_nominal <- function(Ti=NULL,Zr=NULL,Y=NULL,show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.TiZrY_nominal,X=Ti,Y=Zr,Z=Y,f=c(0.01,1,3),
                      labels=c('Ti','Zr','Y'),short=short,
                      show.labels=show.labels,...))
}

xyzplot <- function(json,X=NULL,Y=NULL,Z=NULL,f=rep(1,3),
                    labels=c('X','Y','Z'),short=FALSE,
                    show.labels=FALSE,test.polygons=FALSE,...){
    oldpar <- graphics::par(mar=rep(2,4),xpd=NA)
    ternaryplot(f=f,labels=labels,...)
    if (test.polygons){
        for (pname in names(json$polygons)){
            xyz <-  matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            graphics::polygon(xyz2xy(xyz))
        }
    } else {
        for (lname in names(json$lines)){
            xyz <- matrix(unlist(json$lines[[lname]]),ncol=3,byrow=TRUE)
            graphics::lines(xyz2xy(xyz),lty=lty(json$line_type[[lname]]))
        }
    }
    if (is.null(X) | is.null(Y) | is.null(Z)){
        out <- NULL
    } else {
        XYZ <- sweep(cbind(X,Y,Z),2,f,'*')
        uv <- alr(XYZ,inverse=FALSE)
        XY <- xyz2xy(alr(uv,inverse=TRUE))
        ns <- nrow(XY)
        out <- rep(NA,ns)
        col <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in 1:length(json$polygons)){
            pname <- pnames[i]
            xyz <- matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            matched <- inside(pts=XY,pol=xyz2xy(xyz))
            out[matched] <- json$labels[[pname]]
            col[matched] <- i+1
        }
        ternarypoints(uv,pch=21,bg=col,...)
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xyz <- matrix(unlist(json$label_coords[[lname]]),ncol=3,byrow=TRUE)
            if (short) lab <- lname
            else lab <- json$labels[[lname]]
            graphics::text(xyz2xy(xyz),labels=lab,srt=json$angle[[lname]],pos=1)
        }
    }
    graphics::par(oldpar)
    invisible(out)
}
