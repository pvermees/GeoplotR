# types = nominal type, QDA type, LDA type
XYZhelper <- function(X=NULL,Y=NULL,Z=NULL,type,
                      ternary=TRUE,pch=21,bg=NULL,
                      show.labels=FALSE,short=TRUE,
                      types,nominal,LDA,QDA,xlab,ylab,zlab,
                      f=c(1,1,1),diagram,...){
    if (identical(type[1],types[1])) {
        out <- xyzplot(json=nominal,X=X,Y=Y,Z=Z,f=f,pch=pch,
                       bg=bg,short=short,show.labels=show.labels,...)
    } else {
        uv <- alr(cbind(X,Y,Z))
        quadratic <- identical(type[1],types[2])
        if (quadratic) da <- QDA
        else da <- LDA
        out <- DA(uv=uv,da=da,ternary=ternary,f=f,
                  xlab=xlab,ylab=ylab,zlab=zlab,pch=pch,bg=bg,...)
        plotlabels(diagram=diagram,ternary=ternary,f=f,
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
    invisible(out)
}
XYhelper <- function(X=NULL,Y=NULL,type,
                     ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,
                     short=TRUE,xlim=NULL,ylim=NULL,
                     types,nominal,LDA,QDA,xlab,ylab,zlab,
                     f=c(1,1,1),nf=c(1,1),diagram,m,M,tot=1e6,linear=TRUE,...){
    if (identical(type[1],types[1])){
        good <- !(is.na(X) | is.na(Y))
        if (is.null(xlim)) xlim <- getlimits(x=X[good],m=m[1],M=M[1])*nf[1]
        if (is.null(ylim)) ylim <- getlimits(x=Y[good],m=m[2],M=M[2])*nf[2]
        out <- xyplot(json=nominal,X=nf[1]*X,Y=nf[2]*Y,pch=pch,bg=bg,
                      short=short,show.labels=show.labels,
                      xlim=xlim,ylim=ylim,...)
    } else {
        uv <- alr(cbind(tot-X-Y,X,Y))
        quadratic <- identical(type[1],types[2])
        if (quadratic) da <- QDA
        else da <- LDA
        out <- DA(uv=uv,da=da,D2=TRUE,tot=tot,ternary=ternary,
                  xlab=xlab,ylab=ylab,f=f,pch=pch,bg=bg,xlim=xlim,ylim=ylim,...)
        plotlabels(diagram=diagram,ternary=ternary,f=f,linear=linear,
                   quadratic=quadratic,show.labels=show.labels,short=short)
    }
    invisible(out)
}

#' @title Ti-Zr-Y
#' @description Ti-Zr-Y tectonic discrimination diagram
#' @param Ti vector with Ti concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param Y vector with Y concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Pearce'} for the nominal decision boundaries of Pearce
#'     and Cann (1973).
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with labels for \code{MORB}, \code{IAB} and
#'     \code{OIB}.
#' @references Pearce, J. A., and Cann, J. R., 1973, Tectonic setting
#'     of basic volcanic rocks determined using trace element
#'     analyses: Earth and Planetary Science Letters, v. 19, no. 2,
#'     p. 290-300.
#' @examples
#' data(test,package='GeoplotR')
#' TiZrY(Ti=wtpct2ppm(test[,'TiO2']),
#'       Zr=test[,'Zr'],Y=test[,'Y'],
#'       type='QDA',plot='ternary')
#' @export
TiZrY <- function(Ti=NULL,Zr=NULL,Y=NULL,
                  type=c('Pearce','LDA','QDA'),
                  ternary=TRUE,pch=21,bg=NULL,
                  show.labels=FALSE,short=TRUE,...){
    out <- XYZhelper(X=Ti,Y=Zr,Z=Y,type=type,ternary=ternary,
                     pch=pch,bg=bg,show.labels=show.labels,short=short,
                     types=c('Pearce','QDA','LDA'),
                     nominal=.TiZrY_nominal,LDA=.TiZrY_LDA,QDA=.TiZrY_QDA,
                     xlab='Ti',ylab='Zr',zlab='Y',
                     f=c(1/100,1,3),diagram='TiZrY',...)
    invisible(out)
}

#' @title Nb-Zr-Y
#' @description Nb-Zr-Y tectonic discrimination diagram
#' @param Nb vector with Nb concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param Y vector with Y concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Meschede'} for the nominal decision boundaries of
#'     Meschede (1986).
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with tectonic affinities.
#' @references Meschede, M., 1986. A method of discriminating between
#'     different types of mid-ocean ridge basalts and continental
#'     tholeiites with the Nb-Zr-Y diagram. Chemical geology, 56(3-4),
#'     pp.207-218.
#' @examples
#' data(test,package='GeoplotR')
#' NbZrY(Nb=test[,'Nb'],Zr=test[,'Zr'],Y=test[,'Y'])
#' @export
NbZrY <- function(Nb=NULL,Zr=NULL,Y=NULL,
                  type=c('Meschede','LDA','QDA'),
                  ternary=TRUE,pch=21,bg=NULL,
                  show.labels=FALSE,short=TRUE,...){
    out <- XYZhelper(X=Nb,Y=Zr,Z=Y,type=type,ternary=ternary,
                     pch=pch,bg=bg,show.labels=show.labels,short=short,
                     types=c('Meschede','QDA','LDA'),
                     nominal=.NbZrY_nominal,LDA=.NbZrY_LDA,QDA=.NbZrY_QDA,
                     xlab='Nb',ylab='Zr',zlab='Y',
                     f=c(2,1/4,1),diagram='NbZrY',...)
    invisible(out)
}

#' @title Th-Ta-Hf
#' @description Th-Ta-Hf tectonic discrimination diagram
#' @param Th vector with Th concentrations (ppm)
#' @param Ta vector with Ta concentrations (ppm)
#' @param Hf vector with Hf concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Wood'} for the nominal decision boundaries of Wood
#'     (1980).
#' @param ternary logical. If \code{FALSE}, produces a logratio plot.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with tectonic affinities.
#' @references Wood, D.A., 1980. The application of a Th-Hf-Ta diagram
#'     to problems of tectonomagmatic classification and to
#'     establishing the nature of crustal contamination of basaltic
#'     lavas of the British Tertiary Volcanic Province. Earth and
#'     Planetary Science Letters, 50(1), pp.11-30.
#' @examples
#' data(test,package='GeoplotR')
#' ThTaHf(Th=test[,'Th'],Ta=test[,'Ta'],Hf=test[,'Hf'])
#' @export
ThTaHf <- function(Th=NULL,Ta=NULL,Hf=NULL,
                   type=c('Wood','LDA','QDA'),
                   ternary=TRUE,pch=21,bg=NULL,
                   show.labels=FALSE,short=TRUE,...){
    out <- XYZhelper(X=Hf,Y=Th,Z=Ta,type=type,ternary=ternary,
                     pch=pch,bg=bg,show.labels=show.labels,short=short,
                     types=c('Wood','QDA','LDA'),
                     nominal=.ThTaHf_nominal,LDA=.ThTaHf_LDA,QDA=.ThTaHf_QDA,
                     xlab='Hf',ylab='Th',zlab='Ta',
                     f=c(1,3,3),diagram='ThTaHf',...)
    invisible(out)
}

TiSiSr <- function(Ti=NULL,Si=NULL,Sr=NULL,
                   ternary=TRUE,pch=21,bg=NULL,
                   show.labels=FALSE,short=TRUE,...){
    out <- XYZhelper(X=Ti,Y=Si,Z=Sr,type=type,ternary=ternary,
                     pch=pch,bg=bg,show.labels=show.labels,short=short,
                     types='LDA',LDA=.ThTaHf_LDA,
                     xlab='Ti',ylab='Si',zlab='Sr',
                     f=c(25,1,1000),diagram='TiSiSr',...)
    invisible(out)
}

#' @title Ti-V
#' @description Ti-V tectonic discrimination diagram
#' @param Ti vector with Ti concentrations (ppm)
#' @param V vector with V concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis, or
#'     \code{'Shervais'} for the nominal decision boundaries of
#'     Shervais (1982).
#' @param ternary logical. If \code{TRUE}, produces a ternary diagram.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param show.labels logical. If \code{TRUE}, labels the
#'     discrimination fields.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param short logical. If \code{TRUE}, uses abbreviated labels for
#'     the discrimination fields.
#' @param ... additional arguments for the generic \code{points}
#'     function.
#' @return if \code{type='LDA'} or \code{type='QDA'}, a list with
#'     components \code{class}, \code{posterior} and \code{x};
#'     otherwise a table with labels for \code{MORB}, \code{IAB} and
#'     \code{OIB}.
#' @references Shervais, J.W., 1982. Ti-V plots and the petrogenesis
#'     of modern and ophiolitic lavas. Earth and Planetary Science
#'     Letters, 59(1), pp.101-118.
#' @examples
#' data(test,package='GeoplotR')
#' TiV(Ti=wtpct2ppm(test[,'TiO2']),V=test[,'V'],type='Shervais')
#' @export
TiV <- function(Ti=NULL,V=NULL,type=c('Shervais','LDA','QDA'),
                ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,
                short=TRUE,xlim=NULL,ylim=NULL,...){
    out <- XYhelper(X=Ti,Y=V,type,
                    ternary=ternary,pch=pch,bg=bg,show.labels=show.labels,
                    short=short,xlim=xlim,ylim=ylim,
                    types=c('Shervais','QDA','LDA'),
                    nominal=.TiV_nominal,LDA=.TiV_LDA,QDA=.TiV_QDA,
                    xlab='Ti',ylab='V',f=c(1,100,5000),
                    nf=c(1/1000,1),diagram='TiV',
                    m=c(0,0),M=c(25,600),tot=1e6,linear=TRUE,...)
    invisible(out)
}

#' @title Ti-Zr
#' @description Ti-Zr tectonic discrimination diagram.
#' @param Ti vector with Ti concentrations (ppm)
#' @param Zr vector with Zr concentrations (ppm)
#' @param type either \code{'LDA'} for linear discriminant analysis,
#'     \code{'QDA'} for quadratic discriminant analysis,
#'     \code{'Pearce'} for the nominal decision boundaries of Pearce
#'     and Cann (1973), or \code{'Dilek'} for the nominal decision
#'     boundaries of Dilek and Furnes (2009).
#' @param ternary logical. If \code{TRUE}, produces a ternary diagram.
#'     Only used if \code{type} is \code{LDA} or \code{QDA}.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @return a vector with tectonic affinities
#' @references Pearce, J.A. and Cann, J.R., 1973. Tectonic setting of
#'     basic volcanic rocks determined using trace element
#'     analyses. Earth and Planetary Science Letters, 19(2),
#'     pp. 290-300.
#' 
#' Dilek, Y. and Furnes, H., 2009, Structure and geochemistry of
#' Tethyan ophiolites and their petrogenesis in subduction rollback
#' systems.  Lithos, v. 113, p. 1-20.
#' @examples
#' data(test,package='GeoplotR')
#' ZrTi(Zr=c(50,40,60),Ti=c(1000,4000,6000))
#' @export
ZrTi <- function(Zr=NULL,Ti=NULL,type=c('Pearce','Dilek','LDA','QDA'),
                 ternary=FALSE,pch=21,bg=NULL,show.labels=FALSE,
                 short=TRUE,xlim=NULL,ylim=NULL,...){
    if (type[1]=='Dilek'){
        out <- TiZr_Dilek(Ti=Ti,Zr=Zr,xlim=xlim,ylim=ylim,
                          show.labels=show.labels,short=short,...)
    } else {
        out <- XYhelper(X=Zr,Y=Ti,type,
                        ternary=ternary,pch=pch,bg=bg,show.labels=show.labels,
                        short=short,xlim=xlim,ylim=ylim,
                        types=c('Pearce','QDA','LDA'),
                        nominal=.ZrTi_nominal,LDA=.ZrTi_LDA,QDA=.ZrTi_QDA,
                        xlab='Zr',ylab='Ti',f=c(1,15000,200),
                        nf=c(1,1),diagram='ZrTi',
                        m=c(0,0),M=c(200,15000),tot=1e6,linear=TRUE,...)
    }
    invisible(out)
}
TiZr_Dilek <- function(Ti=NULL,Zr=NULL,xlim=NULL,ylim=NULL,
                       show.labels=TRUE,short=FALSE,...){
    good <- !(is.na(Ti) | is.na(Zr))
    if (is.null(xlim)) xlim <- getlimits(x=Ti[good],m=0,M=10000)
    if (is.null(ylim)) ylim <- getlimits(x=Zr[good],m=0,M=150)
    invisible(xyplot(json=.TiZr,X=Ti,Y=Zr,log='',xlim=xlim,ylim=ylim,
                     show.labels=show.labels,short=short,smooth=TRUE,...))
}

#' @title Vermeesch (2006)
#' @name Vermeesch2006
#' @rdname Vermeesch2006
#' @description discriminant analysis for oceanic basalts
#' @param ternary logical. If \code{TRUE}, produces a ternary
#'     diagram. Otherwise shows the results on a logratio plot.
#' @param pch plot character. See \code{?par} for details.
#' @param bg fill colour for the plot character.
#' @param xlim x-axis limits
#' @param ylim y-axis limits
#' @param show.labels logical. If \code{TRUE}, labels the polygonal
#'     decision fields with text strings.
#' @param short use short labels when using the additional argument
#'     \code{show.labels=TRUE}.
#' @param ... additional arguments to the generic \code{points}
#'     function, may include the logical argument \code{show.labels}
#'     which labels the fields in the diagram.
#' @param Ti vector of Ti-concentrations (in ppm)
#' @param Si vector of Si-concentrations (in ppm)
#' @param Sr vector of Sr-concentrations (in ppm)
#' @return a vector with tectonic affinities
#' @references Vermeesch, 2006. Tectonic discrimination diagrams
#'     revisited: Geochemistry, Geophysics, and Geosystems, 7, Q06017,
#'     doi:10.1029/2005GC001092
#' @examples
#' data(test,package='GeoplotR')
#' TiSiSr(Ti=wtpct2ppm(test$TiO2),Si=wtpct2ppm(test$SiO2),Sr=test$Sr)
#' @export
TiSiSr <- function(Ti=NULL,Si=NULL,Sr=NULL,ternary=FALSE,pch=21,
                   bg=NULL,show.labels=FALSE,short=TRUE,xlim=NULL,
                   ylim=NULL,...){
    uv <- alr(cbind(Ti,Si,Sr))
    f <- c(25,1,1000)
    out <- DA(uv=uv,da=.TiSiSr_LDA,ternary=ternary,f=f,
              xlab='Ti',ylab='Si',zlab='Sr',pch=pch,bg=bg,
              xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='TiSiSr',ternary=ternary,f=f,
               show.labels=show.labels,short=short)
    invisible(out)
}
#' @param Lu vector of Lu-concentrations (in ppm)
#' @param Eu vector of Eu-concentrations (in ppm)
#' @rdname Vermeesch2006
#' @export
LuEuSr <- function(Lu=NULL,Eu=NULL,Sr=NULL,ternary=FALSE,pch=21,
                   bg=NULL,show.labels=FALSE,short=TRUE,xlim=NULL,
                   ylim=NULL,...){
    uv <- alr(cbind(Lu,Eu,Sr))
    f <- c(500,100,1)
    out <- DA(uv=uv,da=.LuEuSr_LDA,ternary=ternary,f=f,
              xlab='Lu',ylab='Eu',zlab='Sr',pch=pch,bg=bg,
              xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='LuEuSr',ternary=ternary,f=f,
               show.labels=show.labels,short=short)
    invisible(out)
}
#' @param V vector of V-concentrations (in ppm)
#' @param Sc vector of Sc-concentrations (in ppm)
#' @rdname Vermeesch2006
#' @export
TiVSc <- function(Ti=NULL,V=NULL,Sc=NULL,ternary=FALSE,pch=21,
                  bg=NULL,show.labels=FALSE,short=TRUE,xlim=NULL,
                  ylim=NULL,...){
    uv <- alr(cbind(Ti,V,Sc))
    f <- c(1,50,250)
    out <- DA(uv=uv,da=.TiVSc_LDA,ternary=ternary,f=f,
              xlab='Ti',ylab='V',zlab='Sc',pch=pch,bg=bg,
              xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='TiVSc',ternary=ternary,f=f,
               show.labels=show.labels,short=short)
    invisible(out)
}
#' @param Nb vector of Nb-concentrations (in ppm)
#' @param Na vector of Na-concentrations (in ppm)
#' @rdname Vermeesch2006
#' @export
NbNaSr <- function(Nb=NULL,Na=NULL,Sr=NULL,ternary=FALSE,pch=21,
                   bg=NULL,show.labels=FALSE,short=TRUE,xlim=NULL,
                   ylim=NULL,...){
    uv <- alr(cbind(Nb,Na,Sr))
    f <- c(2500,1,100)
    out <- DA(uv=uv,da=.NbNaSr_QDA,ternary=ternary,f=f,
              xlab='Nb',ylab='Na',zlab='Sr',pch=pch,bg=bg,
              xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='NbNaSr',ternary=ternary,f=f,
               show.labels=show.labels,short=short)
    invisible(out)
}
#' @param Sm vector of Sm-concentrations (in ppm)
#' @rdname Vermeesch2006
#' @export
TiSmV <- function(Ti=NULL,Sm=NULL,V=NULL,ternary=FALSE,pch=21,
                  bg=NULL,show.labels=FALSE,short=TRUE,xlim=NULL,
                  ylim=NULL,...){
    uv <- alr(cbind(Ti,Sm,V))
    f <- c(1,2500,50)
    out <- DA(uv=uv,da=.TiSmV_QDA,ternary=ternary,f=f,
              xlab='Ti',ylab='Sm',zlab='V',pch=pch,bg=bg,
              xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='TiSmV',ternary=ternary,f=f,
               show.labels=show.labels,short=short)
    invisible(out)
}
#' @param Zr vector of Zr-concentrations (in ppm)
#' @param Y vector of Y-concentrations (in ppm)
#' @rdname Vermeesch2006
#' @export
TiZrYSr <- function(Ti=NULL,Zr=NULL,Y=NULL,Sr=NULL,pch=21,
                    bg=NULL,show.labels=FALSE,short=TRUE,
                    xlim=NULL,ylim=NULL,...){
    uv <- alr(cbind(Ti,Zr,Y,Sr))
    out <- DA_highdim(as.data.frame(uv),da=.TiZrYSr,
                      pch=pch,bg=bg,xlim=xlim,ylim=ylim,...)
    plotlabels(diagram='TiZrYSr',show.labels=show.labels,short=short,raw=TRUE)
    invisible(out)
}
