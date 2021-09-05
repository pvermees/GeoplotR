plotlabels <- function(diagram,ternary=TRUE,f=rep(1,3),
                       quadratic=FALSE,show.labels=TRUE,short=TRUE){
    if (show.labels){
        if (identical(diagram,'AFM')){
            labs <- getlabels(.AFM,short=short)
            if (ternary) uv <- rbind(c(-1.5,-1.5),c(0.5,0.5))
            else uv <- rbind(c(0,-4),c(2,0))
        } else if (identical(diagram,'TiZrY')){
            labs <- getlabels(.TiZrY_nominal,short=short)
            if (quadratic) uv <- rbind(c(-2,-2),c(1,0))
            else uv <- rbind(c(0,-4),c(2,0))
        }
        if (ternary) ternarytext(uv,f=f,labels=labs)
        else graphics::text(x=uv[,1],y=uv[,2],labels=labs)
    }
}

getlabels <- function(json,short=TRUE){
    if (short) return(names(json$labels))
    else return(unlist(json$labels))
}
