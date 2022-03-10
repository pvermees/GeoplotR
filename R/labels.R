plotlabels <- function(diagram,ternary=TRUE,f=rep(1,3),linear=FALSE,
                       quadratic=FALSE,show.labels=TRUE,short=TRUE){
    if (show.labels){
        if (identical(diagram,'AFM')){
            labs <- getlabels(.AFM,short=short)
            if (ternary) uv <- rbind(c(-1.5,-1.5),
                                     c(0.5,0.5))
            else uv <- rbind(c(0,-4),
                             c(2,0))
        } else if (identical(diagram,'TiV')){
            if (quadratic){
                if (short) labs <- c('IAB','IAB','MORB','OIB','OIB')
                else labs <- c('Island Arc','Island Arc','Mid Ocean Ridge',
                               'Ocean Island','Ocean Island')
                uv <- rbind(c(-6,-5),
                            c(-2,-4.5),
                            c(-4.5,-4),
                            c(-4,-7),
                            c(-4,-2.5))
            } else {
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(-5,-5.75),
                            c(-5,-5),
                            c(-5,-7))
            }            
        } else if (identical(diagram,'TiZrY')){
            if (quadratic){
                if (short) labs <- c('IAB','IAB','MORB','OIB','OIB')
                else labs <- c('Island Arc','Island Arc','Mid Ocean Ridge',
                               'Ocean Island','Ocean Island')
                uv <- rbind(c(-6,-5),
                            c(-2,-4.5),
                            c(-4.5,-4),
                            c(-4,-7),
                            c(-4,-2.5))
            } else {
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(-5,-5.75),
                            c(-5,-5),
                            c(-5,-7))
            }
        } else if (identical(diagram,'ZrTi')){
            if (quadratic){
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(-10.5,-6),
                            c(-9.2,-4.6),
                            c(-8,-4))
            } else {
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(-8.5,-6),
                            c(-8.5,-4.5),
                            c(-8.5,-4))
            }            
        }
        if (ternary){
            ternarytext(uv,f=f,labels=labs,xpd=NA)
        } else {
            if (linear) xy <- exp(uv)
            else xy <- uv
            graphics::text(xy[,1:2],labels=labs,xpd=NA)
        }
    }
}

getlabels <- function(json,short=TRUE){
    if (short) return(names(json$labels))
    else return(unlist(json$labels))
}
