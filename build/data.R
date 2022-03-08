setwd('~/Documents/Programming/R/GeoplotR/')
source('R/AFM.R')
source('R/DA.R')
source('R/toolbox.R')

tosave <- NULL

training <- read.csv('inst/training.csv',header=TRUE,check.names=FALSE)
save(training,file="data/training.rda",version=2)

test <- read.csv('inst/test.csv',header=TRUE,check.names=FALSE)
save(test,file="data/test.rda",version=2)

cath <- read.csv('inst/cath.csv',header=TRUE)
save(cath,file="data/cath.rda",version=2)

am <- read.csv('inst/atomicmass.csv',header=FALSE)
.atomicmass <- am[,2]
names(.atomicmass) <- am[,1]
tosave <- c(tosave,'.atomicmass')

.oxides <- read.csv('inst/oxides.csv',header=FALSE,row.names=1)
colnames(.oxides) <- c('cation','ncat','nO')
tosave <- c(tosave,'.oxides')

.TAS <- IsoplotR:::fromJSON(file='inst/TAS.json')
.TAS_volcanic <- IsoplotR:::fromJSON(file='inst/TAS_volcanic.json')
.AnAbOr <- IsoplotR:::fromJSON(file='inst/AnAbOr.json')
.CrY <- IsoplotR:::fromJSON(file='inst/Cr_Y.json')
.LaYb <- IsoplotR:::fromJSON(file='inst/LanYbn_Ybn.json')
.SrY <- IsoplotR:::fromJSON(file='inst/SrY_Y.json')
.ThCo <- IsoplotR:::fromJSON(file='inst/Th_Co.json')
.YNb <- IsoplotR:::fromJSON(file='inst/Pearce_Y-Nb.json')
.YNbRb <- IsoplotR:::fromJSON(file='inst/Pearce_Y+Nb-Rb.json')
.YbTa <- IsoplotR:::fromJSON(file='inst/Pearce_Yb-Ta.json')
.YbTaRb <- IsoplotR:::fromJSON(file='inst/Pearce_Yb+Ta-Rb.json')
.NbLaYb <- IsoplotR:::fromJSON(file='inst/NbLa_LaYb.json')
.ThNbLaYb <- IsoplotR:::fromJSON(file='inst/ThNb_LaYb.json')
.ZrTi <- IsoplotR:::fromJSON(file='inst/Zr_Ti.json')
tosave <- c(tosave,'.TAS','.TAS_volcanic','.AnAbOr','.CrY','.LaYb',
            '.SrY','.ThCo','.YNb','.YNbRb','.YbTa','.YbTaRb',
            '.NbLaYb','.ThNbLaYb','.ZrTi')

.AFM <- IsoplotR:::fromJSON(file='inst/AFM.json')
.QAP <- IsoplotR:::fromJSON(file='inst/QAP.json')
.TiZrY_nominal <- IsoplotR:::fromJSON(file='inst/TiZrY.json')
.TiZrY_LDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
.TiV_LDA <- construct_DA(X='Ti',Y='V',quadratic=FALSE,plot=FALSE)
attributes(.TiZrY_LDA$fit$terms)$.Environment <- NULL
.TiZrY_QDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
.TiV_QDA <- construct_DA(X='Ti',Y='V',quadratic=TRUE,plot=FALSE)
attributes(.TiZrY_QDA$fit$terms)$.Environment <- NULL
tosave <- c(tosave,'.AFM','.QAP','.TiZrY_nominal','.TiZrY_LDA','.TiZrY_QDA')

library(rpart)
treedata_all <- training[c(1,5:55)]
my.control <- rpart.control(xval=10, cp=0, minsplit=1)
unpruned <- rpart(AFFINITY ~ ., data=treedata_all,
                  method="class", control=my.control)
.tectotree_all <- prune(unpruned, cp=0.008)
tosave <- c(tosave,'.tectotree_all')

treedata_HFS <- get_training_data(c("AFFINITY","TiO2","La","Ce","Pr","Nd",
                                    "Sm","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                                    "Lu","Sc","Y","Zr","Nb","Hf","Ta","Pb",
                                    "Th","U","Nd143/Nd144","Sr87/Sr86",
                                    "Pb206/Pb204","Pb207/Pb204","Pb208/Pb204"))
unpruned <- rpart(AFFINITY ~ ., data=treedata_HFS,
                  method="class", control=my.control)
.tectotree_HFS <- prune(unpruned, cp=0.025)
tosave <- c(tosave,'.tectotree_HFS')

num <- c(rep('Ti',23),'Zr','Nb','La','La','Gd',
         'Th','Nb','Th','Th','Nb','Nb','Sr')
den <- c('La','Ce','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
         'Lu', 'Sc','V','Sr','Y','Zr','Nb','Hf','Ta','Th','U','Nb',
         'Th','Sm','Yb','Yb','Ta','La','Yb','U','U','Ta','Zr')
treedata_ratios <- training['AFFINITY']
for (i in 1:length(num)){
    treedata_ratios[paste0(num[i],'/',den[i])] <-
        get_training_data(num[i])/get_training_data(den[i])
}
unpruned <- rpart(AFFINITY ~ ., data=treedata_ratios,
                  method="class", control=my.control)
.tectotree_ratios <- prune(unpruned, cp=0.015)
tosave <- c(tosave,'.tectotree_ratios')

.BF_Fe1 <- list(p=c(9.57,6.98),
                xy0=c(-3.84,2.12))
.BF_Fe2 <- list(p=c(-0.8,-0.3,-1.45,-1,-6,-0.6),
                xyi=c(-1.125,-0.65),
                b=c(-6.0,-0.6),
                d=0.545/sqrt(0.6^2+1))
.BF_Ti <- list(p=c(1,1.4,0,0.65,2.5,-0.35),
               xyi=c(0.5,1.125),
               b=c(2.5,-0.35),
               d=0.55/sqrt(0.35^2+1))
tosave <- c(tosave,'.BF_Fe1','.BF_Fe2','.BF_Ti')

save(list=tosave,file="R/sysdata.rda",version=2)
