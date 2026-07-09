# Self-contained molecular-weight lookup for the oxides used by CIPW.
get_mol_wt <- function(WR){

    mw <- c(
        H=1.008,
        C=12.01,
        O=15.999,
        Na=22.990,
        Mg=24.30,
        Al=26.982,
        Si=28.08,
        P=30.974,
        S=32.0,
        Cl=35.4,
        K=39.098,
        Ca=40.078,
        Ti=47.867,
        Mn=54.938,
        Fe=55.84,
        F=18.998
    )

    molecularWeight <- function(formula){
        ee <- strsplit(paste(formula," ",sep=""),"O")
        n <- as.numeric(gsub("[a-zA-Z]","",ee[[1]]))
        n[is.na(n)] <- 1
        if (length(n) == 1) n[2] <- 0
        atom <- gsub("[0-9 ]","",ee[[1]][1])
        z <- n[1]*mw[atom]+n[2]*mw["O"]
        z <- c(z,n)
        names(z) <- c("MW", "x.atoms", "x.oxygen")
        return(z)
    }

    oxides.simple <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO","MgO",
                       "CaO","Na2O","K2O","H2O","CO2","P2O5","F","S")
    oxides <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO","MgO","CaO",
                "Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")
    
    mol.wt <- sapply(oxides.simple,molecularWeight)
    colnames(mol.wt) <- oxides
    mol.wt
}

millications <- function(x){

    mol.wt <- get_mol_wt()

    x <- data.matrix(x)
    if(ncol(x)==1){
        x <- t(x)
    }
    oxides <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","FeOt","MnO","MgO",
                "CaO","Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")
    oxides <- oxides[oxides%in%colnames(x)]
    ee <- sapply(1:nrow(x),function(f){
        z <- x[f,oxides]/mol.wt["MW",oxides]*mol.wt["x.atoms",oxides]*1000
        return(z)
    })
    milli <- t(ee)
    milli[is.na(milli)] <- 0
    rownames(milli) <- rownames(x)

    # Get ready names of atoms from the oxide names
    .atoms.from.formula <- function(ox){
        z <- gsub("[0-9]","",ox)                          # Remove numbers                        
        z <- sapply((strsplit(z,"O")),paste,collapse="")  # Remove oxygen's "O" and anything beyond
    return(z)
    }

    results <- milli[,oxides,drop=FALSE]
    ee <- paste(.atoms.from.formula(oxides),".m",sep="")
    i <- grep("^Fe.m$",ee)
    if(length(i)==2) ee[i] <- c("Fe3.m","Fe2.m")
    ee <- gsub("^Fet.m$","Fe.m",ee)
    colnames(results) <- ee

    invisible(milli)
}

#' @title CIPW norm
#' @description Cross-Iddings-Pirsson-Washington norm
#' @param wrdata a data frame with whole rock data
#' @param normsum logical. Recast the norm to 100\%?
#' @param cancrinite logical. Is cancrinite present?
#' @param spinel logical. Shall be spinel calculated for rocks with SiO\code{_2 < 45}\%?
#' @param complete.results logical. Remove rows with incomplete results?
#' @return a data frame with the inferred normative mineralogical composition.
#' @references Janousek, V., Farrow, C.M. and Erban, V.,
#' 2006. Interpretation of whole-rock geochemical data in igneous
#' geochemistry: introducing Geochemical Data Toolkit
#' (GCDkit). Journal of Petrology, 47(6), pp.1255-1259.
#'
#' Hutchison, C.S., 1974. Laboratory handbook of petrographic
#' techniques (Vol. 527). New York: Wiley.
#' @author Vojtech Janousek
#' @examples
#' data(cath, package='GeoplotR')
#' norm <- CIPW(cath)
#' @export
CIPW <- function(wrdata,normsum=FALSE,cancrinite=FALSE,
                 spinel=FALSE,complete.results=FALSE){

    CIPW.main <- function(wrdata,normsum=FALSE,cancrinite=FALSE,spinel=FALSE){
        result.names <- c("Q","C","Or","Ab","An","Lc","Ne","Kp","Nc","Ac","Ns","Ks","Di","MgDi",
                          "FeDi","Wo","Hy","En","Fs","Ol","Fo","Fa","Dcs","Mt","Il","Hm","Tn","Pf",
                          "Ru","Ap","Fr","Py","Cc","Sp","MgSp","FeSp","Sum") 
        on.exit(options("show.error.messages"=TRUE))

        ######################################################################
        #                              Main function                         #
        ######################################################################
        CIPWnorm <- function(x){
            names(x) <- c("si","ti","al","fe3","fe2","mn","mg","ca","na","k","H2O","co2","p","fl","s")
            x <- as.list(x)
            si <- x$si
            ti <- x$ti
            al <- x$al
            fe3 <- x$fe3
            fe2 <- x$fe2
            mn <- x$mn
            mg <- x$mg
            ca <- x$ca
            na <- x$na
            k <- x$k
            H2O <- x$H2O
            co2 <- x$co2
            p <- x$p
            fl <- x$fl
            s <- x$s

            fe2 <- fe2 + mn

            norm.names <- c("qtz","cor","ort","ab","an","lc","ne","kp","nc","ac","ns","ks","di",
                            "mgdi","fedi","wo","hy","en","fs","ol","fo","fa","cs","mt","il","hm",
                            "tn","pf","ru","ap","fr","pr","cc","sp","mgsp","fesp","mgr","fer","femg")
            y <- rep(0,length(norm.names))
            names(y) <- norm.names
            y <- as.list(y)

            if (ca >= 10 / 3 * p){
                y$ap <- p
                ca <- ca - y$ap * 10 / 3
            } else {
                y$ap <- 3/10 * ca
                p <- p - y$ap
                ca <- 0
            }

            if (fl >= 2 / 3 * y$ap & !is.na(fl)){
                fl <- fl - 2 / 3 * y$ap
            } else {
                fl <- 0
            }

            if (ca >= .5 * fl){
                y$fr <- .5 * fl
                ca <- ca - y$fr
            } else {
                y$fr <- ca
                fl <- fl - 2 * y$fr
                ca <- 0
            }

            if(!is.na(s)){
                if (fe2 >= .5 * s){
                    y$pr <- .5 * s
                    fe2 <- fe2 - y$pr
                } else {
                    y$pr <- fe2
                    s <- s - 2 * y$pr
                    fe2 <- 0
                }
            }

            if (cancrinite){y$nc <- co2; na <- na - y$nc}

            if(!is.na(co2)){
                if (ca >= co2){
                    y$cc <- co2
                    ca <- ca - y$cc
                    co2 <- 0
                } else {
                    y$cc <- ca
                    co2 <- co2 - y$cc
                    ca <- 0
                }
            }

            if (fe2 >= ti){  
                y$il <- ti
                fe2 <- fe2 - y$il
                ti <- 0
            } else {
                y$il <- fe2
                ti <- ti - y$il
                fe2 <- 0
            }

            if (al >= k){  
                y$ort <- k
                al <- al - y$ort
                si <- si - 6 * y$ort
                k <- 0
            } else {
                y$ort <- al
                k <- k - y$ort
                si <- si - 6 * y$ort
                al <- 0
                y$ks <- k
                si <- si - y$ks
                k <- 0
            }

            if (al >= na){  
                y$ab <- na
                al <- al - y$ab
                si <- si - 6 * y$ab
                na <- 0
            } else {
                y$ab <- al
                na <- na - y$ab
                si <- si - 6 * y$ab
                al <- 0
            }

            if (na >= fe3){
                y$ac <- fe3
                na <- na - y$ac
                fe3 <- 0
                y$ns <- na
                si <- si - 4 * y$ac - y$ns
            } else {
                y$ac <- na
                fe3 <- fe3 - y$ac
                na <- 0
                si <- si - 4 * y$ac
            }

            if (al >= ca){  
                y$an <- ca
                al <- al - y$an
                ca <- 0
                si <- si - 2 * y$an
                y$cor <- al
                al <- 0
            } else {
                y$an <- al
                ca <- ca - y$an
                si <- si - 2 * y$an
                al <- 0
            }

            if (ca >= ti){
                y$tn <- ti
                ca <- ca - y$tn
                si <- si - y$tn
                ti <- 0
            } else {
                y$tn <- ca
                ti <- ti - y$tn
                ca <- 0
                y$ru <- ti
                si <- si - y$tn
                ti <- 0
            }

            if (fe3 >= fe2){
                y$mt <- fe2
                fe3 <- fe3 - y$mt
                fe2 <- 0
                y$hm <- fe3
                fe3 <- 0
            } else {
                y$mt <- fe3
                fe2 <- fe2 - y$mt
                fe3 <- 0
            }

            y$fer <- fe2 / (fe2 + mg)
            y$mgr <- mg / (fe2 + mg)
            y$femg <- fe2 + mg

            if (spinel& si< 45){
                if (y$femg <= y$cor){
                    y$mgsp <- y$mgr * y$femg
                    y$fesp <- y$fer * y$femg
                    y$cor <- y$cor - y$mgsp - y$fesp

                    y$mgsp <- y$mgr * y$cor
                    y$fesp <- y$fer * y$cor
                    y$cor <- 0
                    y$femg <- y$femg - y$mgsp - y$fesp
                }
            }

            if (ca >= y$femg){
                y$di <- y$femg
                ca <- ca - y$femg
                y$wo <- ca
                si <- si - 2 * y$di - y$wo
                ca <- 0
            } else {
                y$di <- ca
                y$femg <- y$femg - ca
                y$hy <- y$femg
                si <- si - 2 * y$di - y$hy
            }

            if (si >= 0){
                y$qtz <- si
                w <- .Ende(y)
                return(w)
            } else {
                y$qtz <- 0
                d <- abs(si)
            }

            if (d <= y$hy / 2){
                y$ol <- d
                y$hy <- y$hy - 2 * d
                w <- .Ende(y)
                return(w)
            } else {
                y$ol <- y$hy / 2
                d <- d - y$hy / 2
                y$hy <- 0
            }

            if (d <= y$tn){
                y$tn <- y$tn - d
                y$pf <- d
                w <- .Ende(y)
                return(w)
            } else {
                y$pf <- y$tn
                d <- d - y$tn
                y$tn <- 0
            }

            if (d <= 4 * y$ab){
                y$ne <- d / 4
                y$ab <- y$ab - d / 4
                w <- .Ende(y)
                return(w)
            } else {
                y$ne <- y$ab
                d <- d - 4 * y$ab
                y$ab <- 0
            }


            if (d <= 2 * y$ort){
                y$lc <- d / 2
                y$ort <- y$ort - d / 2
                w <- .Ende(y)
                return(w)
            } else {
                y$lc <- y$ort
                d <- d - 2 * y$ort
                y$ort <- 0
            }

            if (d < y$wo / 2){
                y$cs <- d
                y$wo <- y$wo - 2 * d
                w <- .Ende(y)
                return(w)
            } else {
                y$cs <- y$wo / 2
                d <- d - y$wo / 2
                y$wo <- 0
            }

            if (d <= y$di){
                y$cs <- y$cs + d / 2
                y$ol <- y$ol + d / 2
                y$di <- y$di - d
                d <- 0
                y$kp <- 0
                w <- .Ende(y)
                return(w)
            } else {
                y$cs <- y$cs + y$di / 2
                y$ol <- y$ol + y$di / 2
                d <- d - y$di
                y$di <- 0
            }

            y$kp <- d / 2
            y$lc <- y$lc - d / 2
            w <- .Ende(y)
            return(w)
        }

        .Ende <- function(y){
            y$en <- y$mgr * y$hy
            y$fs <- y$fer * y$hy
            y$fo <- y$mgr * y$ol
            y$fa <- y$fer * y$ol
            y$mgdi <- y$mgr * y$di
            y$fedi <- y$fer * y$di

            w <- unlist(y)
            w <- w[1:36]

            # Molecular weights calculated
            CIPWweight <- c(60.08480,101.96128,556.66548,524.44902,278.21028,436.49588,284.10982,316.32628,
                            105.98874,462.01034,122.06374,154.28020,1,216.5534,248.09,116.16,1,100.39,
                            131.93,1,140.70,203.78,172.24,231.54,151.75,159.69,196.06,135.98,79.90,
                            336.21,78.08,119.98,100.09,1,142.27,173.81)
            names(CIPWweight) <- names(w)

            w <- w[1:length(CIPWweight)]*CIPWweight

            w[13] <- w[14]+w[15] # di = mgdi + fedi
            w[17] <- w[18]+w[19] # hy = en + fs
            w[20] <- w[21]+w[22] # ol = fo + fa
            w[34] <- w[35]+w[36] # sp = mgsp + fesp

            suma <- sum(w[-c(14,15,18,19,21,22,35,36)])

            # Recast to 100% anhydrous?
            if (normsum){
                w <- w*100/suma
                w[37] <- sum(w[-c(14,15,18,19,21,22,35,36)])
            } else {
                w[37] <- suma
            }
            return(w)
        }

        ######################################################################
        #                              Entry point                           #
        ######################################################################
        oxides <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO","MgO","CaO",
                    "Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")

        # Add missing columns for majors
        ################################
        wrdata <- sapply(oxides,function(f){
            if(!any(colnames(wrdata)==f)){
                return(rep(NA,nrow(wrdata)))
            } else {
                return(wrdata[,f])
            }
        },simplify=TRUE)


        try.it <- function(x){
            CIPWnorm(x)
        }

        if(nrow(wrdata)==0){
            results <- matrix(numeric(0),nrow=0,ncol=length(result.names))
            colnames(results) <- result.names
            return(results)
        }

        MW <- get_mol_wt()[1,]

        results <- sapply(1:nrow(wrdata),function(fff){
            dataset <- wrdata[fff,oxides]/MW[oxides]
            dataset[is.na(dataset)] <- 0
            options(show.error.messages = FALSE)
            res <- try(try.it(dataset))
            options(show.error.messages = TRUE)
            
            if((class(res[1]))!="numeric"){
                cat("Sample",rownames(wrdata)[fff],"- ")
                err <- res[1]
                cat("Error in calculation\n")
                w <- rep(NA,length(result.names))
            } else {
                w <- res
            }
            return(w)
        },simplify=TRUE)

        results <- t(results)
        colnames(results) <- result.names
        rownames(results) <- rownames(wrdata)

        return(results)
    }

    results <- CIPW.main(wrdata,normsum=normsum,cancrinite=cancrinite,spinel=spinel)
    if(!complete.results){
        results <- results[,is.na(match(colnames(results),c("En","Fs","Fo","Fa","MgDi","FeDi")))]
        ee <- apply(results==0,2,sum,na.rm=TRUE)
        results <- results[,ee!=nrow(results)]
        results <- results[!is.na(results[,"Sum"]),]
    }

    invisible(results)
}

#' @title Barth-Niggli's molecular norm
#' @description Calculates the cationic norm of Niggli (1948)
#' @param WR a data frame with whole rock data
#' @return a data frame with the inferred normative mineralogical composition.
#' @references Janousek, V., Farrow, C.M. and Erban, V.,
#' 2006. Interpretation of whole-rock geochemical data in igneous
#' geochemistry: introducing Geochemical Data Toolkit
#' (GCDkit). Journal of Petrology, 47(6), pp.1255-1259.
#'
#' Niggli P (1948) Gesteine und Minerallagerstatten. Birkhauser, Basel, p. 1-540
#' @author Vojtech Janousek
#' @examples
#' data(cath, package='GeoplotR')
#' norm  <-  Catanorm(cath)
#' @export
Catanorm <- function(WR){
    result.names <- c("Q","C","Or","Plag","Ab","An","Lc","Ne","Kp","Ac",
                    "Ns","Ks","Hy","Di","Wo %","En %","Fs %","Ol",
                    "Fo %","Fa %","Cs","Mt","Hm","Il","Tn","Pf",
                    "Ru","Ap","Fr","Py","Cc","Sum")
    on.exit(options("show.error.messages"=TRUE))

    ######################################################################
    #                  Main function
    ######################################################################
    Catanorm.main <- function(x){
        names(x) <- c("si","ti","al","fe3","fe2","mn","mg","ca","na","k","co2","p","fl","s")
        x <- as.list(x)
        si <- x$si
        ti <- x$ti
        al <- x$al
        fe3 <- x$fe3
        fe2 <- x$fe2
        mn <- x$mn
        mg <- x$mg
        ca <- x$ca
        na <- x$na
        k <- x$k
        co2 <- x$co2
        p <- x$p
        fl <- x$fl
        s <- x$s

        norm.names <- c("qtz","cor","ort","plg","ab","an","lc",
                      "ne","kp","ac","ns","ks","hy","di","wo",
                      "en","fs","ol","fo","fa","cs","mt","hm",
                      "il","tn","pf","ru","ap","fr","pr","cc","mgr")
        y <- rep(0,length(norm.names))
        names(y) <- norm.names
        y <- as.list(y)
        
        fe2 <- fe2 + mn
        y$cc <- 2 * co2
        ca <- ca - co2

        if (p < 3 * fl){
            y$ap <- 3 * p; ca <- ca - 1.667 * p; fl <- fl - .33 * p
        }else{
            y$ap <- 2.667 * p + fl; ca <- ca - 1.667 * p; fl <- 0
        }

        y$fr <- 1.5 * fl; ca <- ca - fl / 2; fl <- 0
        y$pr <- 1.5 * s; fe2 <- fe2 - s / 2

        if (ti <= fe2){
            y$il <- 2 * ti; fe2 <- fe2 - ti; ti <- 0
        }else{
            y$il <- 2 * fe2; ti <- ti - fe2; fe2 <- 0
        }

        if (k <= al){
            y$ort <- 5 * k; al <- al - k; si <- si - 3 * k; k <- 0
        }else{
            y$ort <- 5 * al; k <- k - al; si <- si - 3 * al; al <- 0
        }

        y$ks <- 1.5 * k; si <- si - .5 * k; k <- 0

        if (na <= al){
            y$ab <- 5 * na; al <- al - na; si <- si - 3 * na; na <- 0
        }else{
            y$ab <- 5 * al; na <- na - al; si <- si - 3 * al; al <- 0
        }

        if (na <= fe3){
            y$ac <- 4 * na; fe3 <- fe3 - na; si <- si - 2 * na; na <- 0
        }else{
            y$ac <- 4 * fe3; si <- si - 2 * fe3; fe3 <- 0
        }

        y$ns <- 1.5 * na; si <- si - .5 * na; na <- 0

        if (ca <= .5 * al){
            y$an <- 5 * ca; al <- al - 2 * ca; si <- si - 2 * ca; ca <- 0
        }else{
            y$an <- 2.5 * al; ca <- ca - .5 * al; si <- si - al; al <- 0
        }

        if (ti <= ca){
            y$tn <- 3 * ti; ca <- ca - ti; si <- si - ti; ti <- 0
        }else{
            y$tn <- 3 * ca; ti <- ti - ca; si <- si - ca; ca <- 0
        }

        y$ru <- ti; ti <- 0
        y$cor <- al; al <- 0

        if (fe3 <= 2 * fe2){
            y$mt <- 1.5 * fe3; fe2 <- fe2 - fe3 / 2; fe3 <- 0
        }else{
            y$mt <- 3 * fe2; fe3 <- fe3 - 2 * fe2; fe2 <- 0
        }

        y$hm <- fe3; fe3 <- 0
        y$wo <- 2 * ca; si <- si - ca; ca <- 0
        y$en <- 2 * mg; si <- si - mg; mg <- 0
        y$fs <- 2 * fe2; si <- si - fe2; fe2 <- 0
        y$hy <- y$en + y$fs
        y$mgr <- 100 * y$en / (y$en + y$fs)

        if (y$hy < y$wo){
            y$di <- 2 * y$hy; y$wo <- y$wo - y$hy; y$hy <- 0
        }else{
            y$di <- 2 * y$wo; y$hy <- y$hy - y$wo; y$wo <- 0
        }

        if (si >= 0){    y$qtz <- si
            w <- .Ende(y)
            return(w)
        }
        y$qtz <- -si

        if (y$hy >= 4 * y$qtz){
            y$ol <- 3 * y$qtz; y$hy <- y$hy - 4 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$ol <- 3 / 4 * y$hy; y$qtz <- y$qtz - 1 / 4 * y$hy; y$hy <- 0
        }

        if (y$tn >= 3 * y$qtz){
            y$pf <- 2 * y$qtz; y$tn <- y$tn - 3 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$pf <- .667 * y$tn; y$qtz <- y$qtz - 1 / 3 * y$tn; y$tn <- 0
        }

        if (y$ab >= 2.5 * y$qtz){
            y$ne <- 1.5 * y$qtz; y$ab <- y$ab - 2.5 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$ne <- .6 * y$ab; y$qtz <- y$qtz - .4 * y$ab; y$ab <- 0
        }

        if (y$ort >= 5 * y$qtz){
            y$lc <- 4 * y$qtz; y$ort <- y$ort - 5 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$lc <- .8 * y$ort; y$qtz <- y$qtz - .2 * y$ort; y$ort <- 0
        }

        if (y$lc >= 4 * y$qtz){
            y$p <- 3 * y$qtz; y$lc <- y$lc - 4 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$kp <- .75 * y$lc; y$qtz <- y$qtz - .25 * y$lc; y$lc <- 0
        }

        if (y$wo >= 4 * y$qtz){
            y$cs <- 3 * y$qtz; y$wo <- y$wo - 4 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$cs <- .75 * y$wo; y$qtz <- y$qtz - .25 * y$wo; y$wo <- 0
        }

        if (y$di >= 4 * y$qtz){
            y$cs <- y$cs + 1.5 * y$qtz; y$ol <- y$ol + 1.5 * y$qtz; y$di <- y$di - 4 * y$qtz; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }else{
            y$cs <- y$cs + .375 * y$di; y$ol <- y$ol + .375 * y$di; y$qtz <- y$qtz - .25 * y$di; y$di <- 0; y$qtz <- 0
            w <- .Ende(y)
            return(w)
        }
    }

    .Ende <- function(y){
        if(y$di>0){
            y$en <- y$mgr/2
            y$fs <- 50 - y$en
            y$wo <- 50.00
        }else{
            y$en <- 0
            y$fs <- 0
            y$wo <- 0
        }
        if(y$ol>0){
            y$fo <- y$mgr
            y$fa <- 100 - y$fo
        }else{
            y$fo <- 0
            y$fa <- 0
        }

        y$plg <- y$ab+y$an

        w <- unlist(y)
        w <- w[1:31]
        suma <- sum(w[-c(5,6,15,16,17,19,20)])
        w[32] <- suma
        names(w)[32] <- "sum"
        return(w)
    }


    ######################################################################
    #                   Entry point
    ######################################################################
    oxides <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO",
              "MgO","CaO","Na2O","K2O","CO2","P2O5","F","S")

    try.it <- function(x){
        Catanorm.main(x)
    }

    milli <- millications(WR)
    missing.oxides <- oxides[!oxides %in% colnames(milli)]
    if(length(missing.oxides) > 0){
        for(ox in missing.oxides){
            milli <- cbind(milli, stats::setNames(data.frame(rep(0, nrow(milli))), ox))
        }
    }
    milli <- milli[,oxides,drop=FALSE]

    results <- sapply(1:nrow(WR),function(fff){
    
        suma <- sum(milli[fff,oxides],na.rm=TRUE)
        suma <- 100/suma
        wrdata <- milli[fff,oxides]*suma
        wrdata[is.na(wrdata)] <- 0
        options(show.error.messages = FALSE)
        res <- try(try.it(wrdata))
        options(show.error.messages = TRUE)
        
        if((class(res[1]))!="numeric"){
            err <- res[1]
            w <- rep(NA,length(result.names))
        }else{
            w <- res
        }
        return(w)
    },simplify=TRUE)
    
    results <- t(results)
    rownames(results) <- rownames(WR)
    
    results <- results[results[,"sum"]!=0,]
    ee <- apply(results==0,2,sum,na.rm=TRUE)
    results <- results[,ee!=nrow(results)]
        
    if(nrow(results)<nrow(WR)){
        cat("\nNot calculated: \n")
        print(rownames(WR)[is.na(match(rownames(WR),rownames(results)))])
    }

    invisible(results)
}