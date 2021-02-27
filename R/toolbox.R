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
        out <- sweep(log(dat[,-1]),1,log(dat[,1]),'-')
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
rpartpredict <- utils::getFromNamespace("predict.rpart", "rpart")
rparttext <- utils::getFromNamespace("text.rpart", "rpart")

# automatically converts wt% to ppm if necessary
get_training_data <- function(cols){
    out <- NULL
    cnames <- names(training)
    oxides <- c('SiO2','TiO2','Al2O3','Fe2O3','FeO',
                'CaO','MgO','MnO','K2O','Na2O','P2O5')
    elements <- c('Si','Ti','Al','Fe3','Fe2','Ca','Mg','Mn','K','Na','P',
                  'La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho',
                  'Er','Tm','Yb','Lu','Sc','V','Cr','Co','Ni','Cu','Zn',
                  'Ga','Rb','Sr','Y','Zr','Nb','Sn','Cs','Ba','Hf','Ta',
                  'Pb','Th','U')
    isoxide <- (cols %in% oxides)
    iselement <- (cols %in% elements)
    isother <- (cols %in% cnames) & !isoxide & !iselement
    if (any(isoxide)){
        out <- c(out,training[cols[isoxide]])
    }
    if (any(iselement)){
        i <- which(elements %in% cols[iselement])
        if (any(i<12)){ # convert from oxide
            oxide <- oxides[i[i<12]]
            out <- c(out,wtpct2ppm(training[oxide],oxide))
        }
        if (any(i>=12)){
            element <- elements[i[i>=12]]
            out <- c(out,training[element])
        }
    }
    if (any(isother)){
        out <- c(out,training[cols[isother]])
    }
    as.data.frame(out,check.names=FALSE)
}

lty <- function(ltype){
    if (ltype=='dashed') return(2)
    else if (ltype=='dotted') return(3)
    else return(1)
}

# XY = points, xy = polygon
in_polygon <- function(pts,pol){
    if (any(is.na(pts))) return(FALSE)
    ch1 <- chull(pol[,1],pol[,2])
    ch2 <- chull(c(pol[,1],pts[1]),c(pol[,2],pts[2]))
    if (identical(sort(ch1),sort(ch2))) return(TRUE)
    else return(FALSE)
}
