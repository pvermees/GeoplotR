#' @title An-Ab-Or
#' @description An-Ab-Or diagram for tonalite, trondhjemite,
#'     granodiorite, and granite according to Barker (1979).
#' @param An vector with normative anorthite concentrations (wt\%)
#' @param Ab vector with normative albite concentrations (wt\%)
#' @param Or vector with normative orthoclase concentrations (wt\%)
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return a list with labels, indicating whether the input data are
#'     tonalite, trondhjemite, granodiorite, granite.
#' @references Barker, F., 1979, Trondhjemite: Definition,
#'     environment, and hypotheses of origin: p. 1-12, in Barker, F.,
#'     ed., Trondhjemites, Dacites, and Related Rocks, Elsevier,
#'     Amsterdam, 659 p.
#' @examples
#' AnAbOr(An=c(70,75,73),Ab=c(20,10,27),Or=c(10,15,0))
#' @export
AnAbOr <- function(An=NULL,Ab=NULL,Or=NULL,show.labels=TRUE,...){
    invisible(xyzplot(json=.AnAbOr,X=Ab,Y=An,Z=Or,
                      show.labels=show.labels,buffered=TRUE,...))
}

xyzplot <- function(json,X=NULL,Y=NULL,Z=NULL,f=rep(1,3),xyzlab=NULL,
                    show.labels=FALSE,short=FALSE,test.polygons=FALSE,
                    smooth=FALSE,pch=21,bg=NULL,dlwd=1,dcol='black',
                    add=FALSE,neg=FALSE,buffered=FALSE,...){
    oldpar <- graphics::par(mar=rep(2,4),xpd=NA)
    on.exit(graphics::par(oldpar))
    if (is.null(xyzlab)) xyzlab <- json$axis
    if (!add) ternaryplot(f=f,xyzlab=xyzlab,...)
    if (test.polygons){
        pcol <- 2
        for (pname in names(json$polygons)){
            xyz <-  matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            graphics::polygon(xyz2xy(xyz,neg),col=pcol)
            pcol <- pcol+1
        }
    } else {
        for (lname in names(json$lines)){
            xyz <- matrix(unlist(json$lines[[lname]]),ncol=3,byrow=TRUE)
            xy <- xyz2xy(xyz,neg)
            if (smooth) shape <- 1
            else shape <- 0
            graphics::xspline(x=xy[,1],y=xy[,2],shape=shape,border=dcol,
                              lty=lty(json$line_type[[lname]]),lwd=dlwd)
        }
    }
    if (is.null(X) | is.null(Y) | is.null(Z)){
        out <- NULL
    } else {
        mat <- cbind(X,Y,Z)
        XYZ <- sweep(sweep(mat,2,f,'*'),1,rowSums(mat),'/')
        XY <- xyz2xy(XYZ,neg=neg)
        ns <- nrow(XY)
        out <- rep(NA,ns)
        missingbg <- is.null(bg)
        if (missingbg) bg <- rep(1,ns)
        pnames <- names(json$polygons)
        for (i in seq_along(json$polygons)){
            pname <- pnames[i]
            xyz <- matrix(unlist(json$polygons[[pname]]),ncol=3,byrow=TRUE)
            matched <- inside(pts=XY,pol=xyz2xy(xyz,neg),buffered=buffered)
            out[matched] <- ifelse(is.na(out[matched]),json$labels[[pname]],
                                   paste0(out[matched],' + ',json$labels[[pname]]))
            if (missingbg) bg[matched] <- i+1
        }
        points(XY,pch=pch,bg=bg,...)
    }
    if (show.labels){
        for (lname in names(json$labels)){
            xyz <- matrix(unlist(json$label_coords[[lname]]),ncol=3,byrow=TRUE)
            if (short) lab <- lname
            else lab <- json$labels[[lname]]
            a <- angle(json$angle[[lname]])
            graphics::text(xyz2xy(xyz,neg),labels=lab,srt=a,xpd=TRUE)
        }
    }
    invisible(out)
}
