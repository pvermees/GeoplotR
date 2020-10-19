wtpct2ppm <- function(x,oxide){
    oxides <- c('SiO2','TiO2','Al2O3','Fe2O3','FeO',
                'CaO','MgO','MnO','K2O','Na2O','P2O5')
    cations <- c('Si','Ti','Al','Fe3','Fe2',
                 'Ca','Mg','Mn','K','Na','P')
    ncat <- c(1,1,2,2,1,1,1,1,2,2,2)
    nO <- c(2,2,3,3,1,1,1,1,1,1,5)
    i <- which(oxides %in% oxide)
    num <- .atomicmass[cations[i]]
    den <- ncat[i]*.atomicmass[cations[i]] + nO[i]*.atomicmass['O']
    1e4*x*num/den
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
