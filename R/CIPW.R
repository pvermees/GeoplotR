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
#' data(test, package='GeoplotR')
#' norm <- CIPW(test)
#' @export
CIPW <- function(wrdata,normsum=FALSE,cancrinite=FALSE,
                 spinel=FALSE,complete.results=FALSE){
    # Self-contained molecular-weight lookup for the oxides used by CIPW.
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
        if(length(n)==1)n[2] <- 0
        atom <- gsub("[0-9 ]","",ee[[1]][1])
        z <- n[1]*mw[atom]+n[2]*mw["O"]
        z <- c(z,n)
        names(z) <- c("MW","x.atoms","x.oxygen")
        return(z)
    }

    oxides.simple <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO","MgO",
                       "CaO","Na2O","K2O","H2O","CO2","P2O5","F","S")
    oxides <- c("SiO2","TiO2","Al2O3","Fe2O3","FeO","MnO","MgO","CaO",
                "Na2O","K2O","H2O.PLUS","CO2","P2O5","F","S")
    mol.wt <- sapply(oxides.simple,molecularWeight)
    colnames(mol.wt) <- oxides
    MW <- mol.wt[1,]

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
            attach(x,warn.conflicts=FALSE)
            on.exit(detach("x"), add=TRUE)

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
                                        #cat(err)
                w <- rep(NA,length(result.names))
                                        #w <- NA
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

                  
