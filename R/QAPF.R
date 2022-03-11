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
#' @examples
#' QAP(Q=c(10,20),A=c(10,30),P=c(80,50))
#' @return a vector with rock names
#' @export
QAPF <- function(Q=NULL,A=NULL,P=NULL,F=NULL){
    
}
#' @rdname QAPF
#' @export
QAP <- function(Q=NULL,A=NULL,P=NULL,pch=21,bg=NULL,
                show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.QAP,X=Q,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}
#' @rdname QAPF
#' @export
FAP <- function(F=NULL,A=NULL,P=NULL,pch=21,bg=NULL,
                show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.FAP,X=F,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}
