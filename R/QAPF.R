QAPF <- function(Q=NULL,A=NULL,P=NULL,F=NULL){
    
}

#' @title QAP diagram
#' @description Streckeisen (1967)'s Quartz - Alkali feldspar -
#'     Plagioclase feldspar diagram for classification of felsic and
#'     intermediate plutonic rocks.
#' @param Q vector with quartz concentrations (\%)
#' @param A vector with alkali feldspar concentrations (\%)
#' @param P vector with plagioclase feldspar concentrations (\%)
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour of the plot symbols.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with rock types
#' @references Streckeisen, A., 1967, Classification and nomenclature
#'     of igneous rocks.  Nues Jarbuch fur Mineralogie Abhandlungen,
#'     v. 107, p. 144-240.
#' @examples
#' QAP(Q=c(40,30),A=c(20,40),P=c(40,30))
#' @export
QAP <- function(Q=NULL,A=NULL,P=NULL,pch=21,bg=NULL,
                show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.QAP,X=Q,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}

#' @title QAP diagram
#' @description Streckeisen (1967)'s Quartz - Alkali feldspar -
#'     Plagioclase feldspar diagram for classification of felsic and
#'     intermediate plutonic rocks.
#' @param F vector with feldspathoid concentrations (\%)
#' @param A vector with alkali feldspar concentrations (\%)
#' @param P vector with plagioclase feldspar concentrations (\%)
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour of the plot symbols.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields on the plot.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with rock types
#' @references Streckeisen, A., 1967, Classification and nomenclature
#'     of igneous rocks.  Nues Jarbuch fur Mineralogie Abhandlungen,
#'     v. 107, p. 144-240.
#' @examples
#' FAP(F=c(40,30),A=c(20,40),P=c(40,30))
#' @export
FAP <- function(F=NULL,A=NULL,P=NULL,pch=21,bg=NULL,
                show.labels=TRUE,short=TRUE,...){
    invisible(xyzplot(json=.FAP,X=F,Y=A,Z=P,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,...))
}
