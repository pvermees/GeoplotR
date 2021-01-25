#' tectonic discrimination by decision tree
#' @description Classification And Regression Trees (CART) of basalts
#'     from ocean islands (OIB), mid-oceanic ridges (MORB) and island
#'     arcs (IAB).
#' @param dat a data frame or matrix with the following columns
#'     (missing data are allowed):
#'
#' If \code{option=1}: \code{SiO2} (as wt\%) or \code{Si} (as ppm),
#' \code{TiO2} (as wt\%) or \code{Ti} (as ppm), \code{CaO} (as wt\%)
#' or \code{Ca} (as ppm), \code{Fe2O3} (as wt\%) or \code{Fe3} (as
#' ppm), \code{MgO} (as wt\%) or \code{Mg} (as ppm), \code{K2O} (as
#' wt\%) or \code{K} as ppm, \code{La}, \code{Pr}, \code{Nd},
#' \code{Sm}, \code{Gd}, \code{Tb}, \code{Yb}, \code{Lu}, \code{V},
#' \code{Ni}, \code{Rb}, \code{Sr}, \code{Y}, \code{Hf}, \code{Th}
#' (all as ppm), \code{Sr87/Sr86} and \code{Pb206/Pb204}.
#'
#' If \code{option=2}: \code{TiO2} (as wt\%) or \code{Ti} (as ppm),
#' \code{La}, \code{Ce}, \code{Pr}, \code{Nd}, \code{Sm}, \code{Gd},
#' \code{Tb}, \code{Dy}, \code{Ho}, \code{Er}, \code{Tm}, \code{Yb},
#' \code{Lu}, \code{Sc}, \code{Y}, \code{Zr}, \code{Nb}, \code{Hf},
#' \code{Ta}, \code{Pb}, \code{Th}, \code{U} (all as ppm),
#' {Nd143/Nd144}, \code{Sr87/Sr86}, \code{Pb206/Pb204},
#' \code{Pb207/Pb204} and \code{Pb208/Pb204}.
#'
#' If \code{option=3}: \code{TiO2} (as wt\%) or \code{Ti} (as ppm),
#' \code{La}, \code{Sm}, \code{Nd}, \code{Gd}, \code{Yb}, \code{Sc},
#' \code{V}, \code{Sr}, {Y}, \code{Zr}, \code{Nb}, \code{Th}, \code{U}
#' (all as ppm).
#' @param option numerical. If \code{option=1}, uses all major and
#'     trace element concentrations and isotopic ratios, if
#'     \code{option=2}, uses high field strength element
#'     concentrations and isotopic ratios, if \code{option=3}, uses
#'     high field strength element ratios.
#' @export
cart <- function(dat,option=1){
    cdat <- dat4cart(dat,option=option)
    if (option==1) tree <- .tectotree_all
    else if (option==2) tree <- .tectotree_HFS
    else if (option==3) tree <- .tectotree_ratios
    else stop('Illegal option provided to cart function.')
    out <- rpart:::predict.rpart(tree,newdata=cdat)
    invisible(out)
}

dat4cart <- function(dat,option=1){
    oxides <- c('SiO2','TiO2','Al2O3','Fe2O3','FeO',
                'CaO','MgO','MnO','K2O','Na2O','P2O5')
    elements <- c('Si','Ti','Al','Fe3','Fe2','Ca','Mg','Mn','K','Na','P')
    dnames <- names(dat)
    out <- list()
    for (dname in dnames){
        if (dname%in%oxides & option==3){
            i <- which(oxides%in%dname)
            oxide <- oxides[i]
            element <- elements[i]
            out[[element]] <- wtpct2ppm(dat[,dname],oxide)
        } else if (dname%in%elements & option<3){
            i <- which(elements%in%dname)
            oxide <- oxides[i]
            element <- elements[i]
            out[[oxide]] <- ppm2wtpct(dat[,dname],oxide)
        } else {
            out[[dname]] <- dat[,dname]
        }
    }
    if (option==3){
        num <- c(rep('Ti',23),'Zr','Nb','La','La','Gd',
                 'Th','Nb','Th','Th','Nb','Nb','Sr')
        den <- c('La','Ce','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
                 'Lu', 'Sc','V','Sr','Y','Zr','Nb','Hf','Ta','Th','U','Nb',
                 'Th','Sm','Yb','Yb','Ta','La','Yb','U','U','Ta','Zr')
        ratios <- list()
        for (i in 1:length(num)){
            if (num[i]%in%names(out) & den[i]%in%names(out)){
                ratios[[paste0(num[i],'/',den[i])]] <-
                    out[[num[i]]]/out[[den[i]]]
            }
        }
        out <- ratios
    }
    as.data.frame(out,check.names=FALSE)
}

