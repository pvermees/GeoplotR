#' @title concentration unit conversion
#' @description Converts from percentage oxide by weight to parts per
#'     million of the element and vice versa. \code{wtpct2ppm} and
#'     \code{ppm2wtpct} are special cases.
#' @param x scalar or vector of wt\% or ppm values
#' @param oxide Oxide of the element; one of \code{'SiO2'},
#'     \code{'TiO2'}, \code{'Al2O3'}, \code{'Fe2O3'}, \code{'FeO'},
#'     \code{'CaO'}, \code{'MgO'}, \code{'MnO'}, \code{'K2O'},
#'     \code{'Na2O'}, \code{'P2O5'}, \code{'Y2O3'}, \code{'ZrO2'}
#' @param wtpct2ppm logical. If \code{TRUE}, converts wt\% to ppm,
#'     otherwise converts ppm to wt\%.
#' @return converted values
#' @examples
#' ppm <- wtpct2ppm(0.308,'TiO2')
#' wtpct <- ppm2wtpct(ppm,'TiO2')
#' @export
conconv <- function(x,oxide,wtpct2ppm=TRUE){
    cation <- as.character(.oxides[oxide,'cation'])
    ncat <- .oxides[oxide,'ncat']
    nO <- .oxides[oxide,'nO']
    num <- .atomicmass[cation]
    den <- ncat*.atomicmass[cation] + nO*.atomicmass['O']
    if (wtpct2ppm) out <- 1e4*x*num/den
    else out <- x*1e-4*den/num
    out
}
#' @rdname conconv
#' @examples
#' ppm <- wtpct2ppm(0.308,'TiO2')
#' @export
wtpct2ppm <- function(x,oxide){
    conconv(x=x,oxide=oxide,wtpct2ppm=TRUE)
}
#' @rdname conconv
#' @examples
#' wtpct <- ppm2wtpct(ppm,'TiO2')
#' @export
ppm2wtpct <- function(x,oxide){
    conconv(x=x,oxide=oxide,wtpct2ppm=FALSE)
}

alr <- function(dat,inverse=FALSE){
    if (inverse){
        num <- cbind(1,exp(dat))
        den <- 1+rowSums(exp(dat),na.rm=TRUE)
        out <- num/den        
    } else {
        out <- log(dat[,-1])-log(dat[,1])
    }
    out
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
