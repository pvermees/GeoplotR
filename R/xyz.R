AlFeTiMg <- function(Al,Fe,Ti,Mg,...){
    xyplot(x=Al+Fe,y=Ti+Mg,coords=.AlFeTiMg$coords,...)
}

AnAbOr <- function(An,Ab,Or,...){
    ternaryplot(labels=c('An','Ab','Or'))
    lnames <- names(.AnAbOr$polygons)
    coords <- .AnAbOr$polygons
    for (lname in lnames){
#        if (.AnAbOr$type[[lname]]=="dashed") lty <- 2
#        else if (.AnAbOr$type[[lname]]=="dotted") lty <- 3
#        else lty <- 1
        xyz <- matrix(unlist(coords[[lname]]),ncol=3,byrow=TRUE)
        graphics::lines(xyz2xy(xyz),lty=1)
    }
    lnames <- names(.AnAbOr$labels)
    coords <- .AnAbOr$labels
    angles <- .AnAbOr$angle
    for (lname in lnames){
        graphics::text(xyz2xy(coords[[lname]]),
                       labels=lname,srt=angles[[lname]])
    }
}
