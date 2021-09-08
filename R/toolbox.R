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
        out <- sweep(num,1,den,'/')
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

# This uses the ray-casting algorithm to decide whether the point is inside
# the given polygon. See https://en.wikipedia.org/wiki/Point_in_polygon
inside <- function(pts,pol,log=''){
    nv <- nrow(pol)
    if (identical(pol[1,],pol[nv,])){
        pol <- pol[-1,] # remove the duplicate vertex
        nv <- nv - 1
    }
    if (class(pts)=='matrix'){
        np <- nrow(pts)
        x <- pts[,1]
        y <- pts[,2]
    } else {
        np <- 1
        x <- pts[1]
        y <- pts[2]
    }
    if (log%in%c('x','xy')){
        pol[,1] <- log(pol[,1])
        x <- log(x)
    }
    if (log%in%c('x','xy')){
        pol[,2] <- log(pol[,2])
        y <- log(y)
    }
    igood <- which(!(is.na(x)|is.na(y)))
    out <- rep(FALSE,np)
    for (i in 1:nv){
        j <- i %% nv + 1
        xp0 <- pol[i,1]
        yp0 <- pol[i,2]
        xp1 <- pol[j,1]
        yp1 <- pol[j,2]
        crosses <- (yp0 <= y) & (yp1 > y) | (yp1 <= y) & (yp0 > y)
        if (any(crosses[igood])){
            icrosses <- igood[which(crosses[igood])]
            cross <- (xp1 - xp0) * (y[icrosses] - yp0) / (yp1 - yp0) + xp0
            change <- icrosses[cross < x[icrosses]]
            out[change] <- !out[change]
        }
    }
    out
}

angle <- function(a){
    usr <- graphics::par('usr')
    dx <- usr[2]-usr[1]
    dy <- usr[4]-usr[3]
    atan(tan(a*pi/180)*dx/dy)*180/pi
}

getlimits <- function(x,m,M){
    if (is.null(x)){
        out <- c(m,M)
    } else if (m<M){
        out <- range(x,na.rm=TRUE)
        out[1] <- min(out[1],m)
        out[2] <- max(out[2],M)
    } else {
        out <- range(x,na.rm=TRUE)
        out[1] <- max(out[2],m)
        out[2] <- min(out[1],M)        
    }
    invisible(out)
}

# helper function to generate .json files
getXYZ <- function(x,yz){
    X <- x
    Z <- (1-x)/(1+yz)
    Y <- 1-X-Z
    100*c(X,Y,Z)
}
