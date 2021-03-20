xytest <- function(fname,xlim=NA,ylim=NA,log='',short=FALSE,polygons=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    if (any(is.na(xlim)) & any(is.na(ylim))){
        xy <- NULL
        for (pname in names(json$polygons)){
            XY <- matrix(unlist(json$polygons[[pname]]),ncol=2,byrow=TRUE)
            xy <-  rbind(xy,XY)
        }
        xlim <- range(xy[,1])
        ylim <- range(xy[,2])
    }
    out <- xyplot(json,xlim=xlim,ylim=ylim,show.labels=TRUE,
                  short=short,test.polygons=polygons,log=log,...)
    invisible(out)
}

xyztest <- function(fname,short=FALSE,polygons=FALSE,...){
    json <- IsoplotR:::fromJSON(file=fname)
    out <- xyzplot(json,short=short,show.labels=TRUE,
                   test.polygons=polygons,...)
    invisible(out)
}
