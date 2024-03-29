plotlabels <- function(diagram,ternary=TRUE,f=rep(1,3),linear=FALSE,
                       quadratic=FALSE,show.labels=TRUE,short=TRUE,
                       raw=FALSE){
    if (show.labels){
        if (identical(diagram,'AFM')){
            labs <- getlabels(.AFM,short=short)
            if (ternary) uv <- rbind(c(-1.5,-1.5),
                                     c(0.5,0.5))
            else uv <- rbind(c(0,-4),
                             c(2,0))
        } else if (identical(diagram,'TiV')){
            if (short) labs <- c('IAB','MORB','OIB','OIB','IAB')
            else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
            uv <- rbind(c(8.0,5.5),
                        c(9.1,5.5),
                        c(9.9,5.5),
                        c(9.6,8.5),
                        c(12.8,5.9))
            if (ternary) uv <- uv - 14 # from ppm to fraction
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
        } else if (identical(diagram,'NbZrY')){
            if (quadratic){
                if (short) labs <- c('IAB','IAB','MORB','OIB')
                else labs <- c('Island Arc','Island Arc',
                               'Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(1,1),
                            c(3.5,1),
                            c(3,2),
                            c(2,0))
            } else {
                if (short) labs <- c('IAB','OIB','MORB')
                else labs <- c('Island Arc','Ocean Island','Mid Ocean Ridge')
                uv <- rbind(c(3.5,2),
                            c(2,0),
                            c(2,1.5))
            }
        } else if (identical(diagram,'ThTaHf')){
            if (quadratic){
                if (short) labs <- c('IAB','IAB','MORB','OIB','OIB')
                else labs <- c('Island Arc','Island Arc',
                               'Mid Ocean Ridge','Ocean Island','Ocean Island')
                uv <- rbind(c(-1,-2.25),
                            c(-3,-1),
                            c(-1,-1),
                            c(-2.5,-2.25),
                            c(1.55,1.3))
            } else {
                if (short) labs <- c('IAB','OIB','MORB')
                else labs <- c('Island Arc','Ocean Island','Mid Ocean Ridge')
                uv <- rbind(c(-1,-2.5),
                            c(-1,-1),
                            c(-2.5,-2.25))
            }
        } else if (identical(diagram,'TiSiSr')){
            if (short) labs <- c('IAB','OIB','MORB')
            else labs <- c('Island Arc','Ocean Island','Mid Ocean Ridge')
            uv <- rbind(c(3.75,-3),
                        c(2.5,-3.5),
                        c(3.25,-5))
        } else if (identical(diagram,'LuEuSr')){
            if (short) labs <- c('IAB','OIB','MORB')
            else labs <- c('Island Arc','Ocean Island','Mid Ocean Ridge')
            uv <- rbind(c(1,7),
                        c(2,7),
                        c(1,5))
        } else if (identical(diagram,'TiVSc')){
            if (short) labs <- c('IAB','OIB','MORB')
            else labs <- c('Island Arc','Ocean Island','Mid Ocean Ridge')
            uv <- rbind(c(-2.5,-5.5),
                        c(-4,-6.5),
                        c(-4,-5.5))
        } else if (identical(diagram,'NbNaSr')){
            if (short) labs <- c('IAB','IAB','OIB','MORB')
            else labs <- c('Island Arc','Island Arc',
                           'Ocean Island','Mid Ocean Ridge')
            uv <- rbind(c(8,4.5),
                        c(7,1),
                        c(5.5,2.5),
                        c(7,2.5))
        } else if (identical(diagram,'TiSmV')){
            if (short) labs <- c('IAB','IAB','OIB','MORB')
            else labs <- c('Island Arc','Island Arc',
                           'Ocean Island','Mid Ocean Ridge')
            uv <- rbind(c(-9.5,-6),
                        c(-7,-2.5),
                        c(-7.75,-4.5),
                        c(-7.9,-3.25))
        } else if (identical(diagram,'ZrTi')){
            if (quadratic){
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(3.5,8.5),
                            c(4.8,9.2),
                            c(5.2,9.9))
            } else {
                if (short) labs <- c('IAB','MORB','OIB')
                else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
                uv <- rbind(c(5.5,8.5),
                            c(5.5,9.5),
                            c(5.5,9.9))
            }
            if (ternary) uv <- uv - 14 # from ppm to fraction
        } else if (identical(diagram,'TiZrYSr')){
            if (short) labs <- c('IAB','MORB','OIB')
            else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
            uv <- rbind(c(0,2),
                        c(-2,-1),
                        c(1,-1))
        } else if (identical(diagram,'Pearce1976')){
            if (short) labs <- c('IAB','MORB','OIB')
            else labs <- c('Island Arc','Mid Ocean Ridge','Ocean Island')
            uv <- rbind(c(-1,2.5),
                        c(0,-2),
                        c(1,2))
        }
        if (raw){
            graphics::text(uv,labels=labs,xpd=NA)
        } else if (ternary){
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
