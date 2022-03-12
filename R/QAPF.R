#' @title QAPF diagram
#' @name QAPF
#' @rdname QAPF
#' @description Streckeisen (1967)'s Quartz - Alkali feldspar -
#'     Plagioclase feldspar diagram for classification of felsic and
#'     intermediate plutonic rocks.
#' @param Q quartz (\%)
#' @param A alkali feldspar (\%)
#' @param P plagioclase feldspar (\%)
#' @param F feldspathoid (\%)
#' @param volcanic logical. If \code{TRUE}, uses the volcanic
#'     classification of Le Maitre et al. (2004).
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour of the plot symbols.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @references Streckeisen, A., 1967, Classification and nomenclature
#'     of igneous rocks.  Nues Jarbuch fur Mineralogie Abhandlungen,
#'     v. 107, p. 144-240.
#'
#' Le Maitre, R. W., Streckeisen, A., Zanettin, B., Bas, M. J. L.,
#' Bonin, B., and Bateman, P., 2004, Igneous Rocks: A Classification
#' and Glossary of Terms: Cambridge University Press.
#' @examples
#'    out <- QAPF(Q=c(50,40,0,0,80,60),
#'                A=c(10,20,50,0,0,20),
#'                P=c(40,40,50,20,20,20),
#'                F=c(0,0,0,80,0,0))
#' @return a vector with rock names
#' @export
QAPF <- function(Q=NULL,A=NULL,P=NULL,F=NULL,volcanic=FALSE,
                 pch=21,bg=NULL,show.labels=TRUE,short=TRUE,...){
    diamondplot(labels=c('A','Q','P','F'))
    if (volcanic) lines(c(0,0.65),c(0,0))
    else lines(c(0,1),c(0,0))
    out <- rep(NA,max(length(Q),length(F)))
    top <- (Q>=0 & F==0)
    if (volcanic){
        topjson <- .QAP_volcanic
        bottomjson <- .FAP_volcanic
        # the following 6 lines merge the BA fields of QAP and FAP
        bottomjson$labels[['BA']] <- NULL
        FAPBA <- (F>0 & P>65 & F<10)
        out[FAPBA] <- "basalt/andesite"
        col <- which(topjson$labels %in% "basalt/andesite")+1
        points(xyz2xy(cbind(F,A,P)[FAPBA,],neg=TRUE),pch=pch,bg=col,...)
        bottom <- (F>0 & !FAPBA)
    } else {
        topjson <- .QAP
        bottomjson <- .FAP
        bottom <- (F>0)
    }
    out[top] <- xyzplot(json=topjson,X=Q[top],Y=A[top],Z=P[top],
                        pch=pch,bg=bg,short=short,
                        show.labels=show.labels,add=TRUE,
                        buffered=TRUE,...)
    out[bottom] <- xyzplot(json=bottomjson,X=F[bottom],Y=A[bottom],Z=P[bottom],
                           pch=pch,bg=bg,short=short,show.labels=show.labels,
                           add=TRUE,neg=TRUE,buffered=TRUE,...)
    invisible(out)
}
#' @rdname QAPF
#' @export
QAP <- function(Q=NULL,A=NULL,P=NULL,volcanic=FALSE,
                pch=21,bg=NULL,show.labels=TRUE,short=TRUE,...){
    if (volcanic) son <- .QAP_volcanic
    else json <- .QAP
    invisible(xyzplot(json=json,X=Q,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}
#' @rdname QAPF
#' @export
FAP <- function(F=NULL,A=NULL,P=NULL,pch=21,bg=NULL,
                show.labels=TRUE,short=TRUE,...){
    if (volcanic) son <- .FAP_volcanic
    else json <- .FAP
    invisible(xyzplot(json=json,X=F,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}
