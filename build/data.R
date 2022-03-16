setwd('~/Documents/Programming/R/GeoplotR/')
source('R/AFM.R')
source('R/DA.R')
source('R/toolbox.R')

tosave <- NULL

message('load csv files')

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

message('read json files')

.TAS <- IsoplotR:::fromJSON(file='inst/TAS.json')
.TAS_plutonic <- IsoplotR:::fromJSON(file='inst/TAS_plutonic.json')
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
.TiZr <- IsoplotR:::fromJSON(file='inst/TiZr.json')
.AFM <- IsoplotR:::fromJSON(file='inst/AFM.json')
.QAP <- IsoplotR:::fromJSON(file='inst/QAP.json')
.FAP <- IsoplotR:::fromJSON(file='inst/FAP.json')
.QAP_volcanic <- IsoplotR:::fromJSON(file='inst/QAP_volcanic.json')
.FAP_volcanic <- IsoplotR:::fromJSON(file='inst/FAP_volcanic.json')
tosave <- c(tosave,'.TAS','.TAS_plutonic','.AnAbOr','.CrY','.LaYb',
            '.SrY','.ThCo','.YNb','.YNbRb','.YbTa','.YbTaRb',
            '.NbLaYb','.ThNbLaYb','.TiZr','.AFM',
            '.QAP','.FAP','.QAP_volcanic','.FAP_volcanic')

message('Build TiZrY DA')
.TiZrY_nominal <- IsoplotR:::fromJSON(file='inst/TiZrY.json')
.TiZrY_LDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
.TiZrY_QDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
attributes(.TiZrY_LDA$fit$terms)$.Environment <- NULL
attributes(.TiZrY_QDA$fit$terms)$.Environment <- NULL
message('Build NbZrY DA')
.NbZrY_nominal <- IsoplotR:::fromJSON(file='inst/NbZrY.json')
.NbZrY_LDA <- construct_DA(X='Nb',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
.NbZrY_QDA <- construct_DA(X='Nb',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
attributes(.NbZrY_LDA$fit$terms)$.Environment <- NULL
attributes(.NbZrY_QDA$fit$terms)$.Environment <- NULL
message('Build ThTaHf DA')
.ThTaHf_nominal <- IsoplotR:::fromJSON(file='inst/ThTaHf.json')
.ThTaHf_LDA <- construct_DA(X='Hf',Y='Th',Z='Ta',quadratic=FALSE,plot=FALSE)
.ThTaHf_QDA <- construct_DA(X='Hf',Y='Th',Z='Ta',quadratic=TRUE,plot=FALSE)
attributes(.ThTaHf_LDA$fit$terms)$.Environment <- NULL
attributes(.ThTaHf_QDA$fit$terms)$.Environment <- NULL
message('Build TiSiSr DA')
.TiSiSr_LDA <- construct_DA(X='Ti',Y='Si',Z='Sr',quadratic=FALSE,plot=FALSE)
attributes(.TiSiSr_LDA$fit$terms)$.Environment <- NULL
message('Build LuEuSr DA')
.LuEuSr_LDA <- construct_DA(X='Lu',Y='Eu',Z='Sr',quadratic=FALSE,plot=FALSE)
attributes(.LuEuSr_LDA$fit$terms)$.Environment <- NULL
message('Build TiVSc DA')
.TiVSc_LDA <- construct_DA(X='Ti',Y='V',Z='Sc',quadratic=FALSE,plot=FALSE)
attributes(.TiVSc_LDA$fit$terms)$.Environment <- NULL
message('Build NbNaSr DA')
.NbNaSr_QDA <- construct_DA(X='Nb',Y='Na',Z='Sr',quadratic=TRUE,plot=FALSE)
attributes(.NbNaSr_QDA$fit$terms)$.Environment <- NULL
message('Build TiSmV DA')
.TiSmV_QDA <- construct_DA(X='Ti',Y='Sm',Z='V',quadratic=TRUE,plot=FALSE)
attributes(.TiSmV_QDA$fit$terms)$.Environment <- NULL
message('Build TiV DA')
.TiV_nominal <- IsoplotR:::fromJSON(file='inst/TiV.json')
.TiV_LDA <- construct_DA(X='Ti',Y='V',quadratic=FALSE,plot=FALSE)
.TiV_QDA <- construct_DA(X='Ti',Y='V',quadratic=TRUE,plot=FALSE)
attributes(.TiV_LDA$fit$terms)$.Environment <- NULL
attributes(.TiV_QDA$fit$terms)$.Environment <- NULL
message('Build ZrTi DA')
.ZrTi_nominal <- IsoplotR:::fromJSON(file='inst/ZrTi.json')
.ZrTi_LDA <- construct_DA(X='Zr',Y='Ti',quadratic=FALSE,plot=FALSE)
.ZrTi_QDA <- construct_DA(X='Zr',Y='Ti',quadratic=TRUE,plot=FALSE)
attributes(.ZrTi_LDA$fit$terms)$.Environment <- NULL
attributes(.ZrTi_QDA$fit$terms)$.Environment <- NULL
message('Build TiZrYSr DA')
.TiZrYSr <- construct_DA_highdim(c('Ti','Zr','Y','Sr'),quadratic=FALSE,plot=FALSE)
attributes(.TiZrYSr$fit$terms)$.Environment <- NULL
tosave <- c(tosave,
            '.TiSiSr_LDA','.LuEuSr_LDA','.TiVSc_LDA',
            '.NbNaSr_QDA','.TiSmV_QDA',
            '.ThTaHf_nominal','.ThTaHf_LDA','.ThTaHf_QDA',
            '.TiZrY_nominal','.TiZrY_LDA','.TiZrY_QDA',
            '.NbZrY_nominal','.NbZrY_LDA','.NbZrY_QDA',
            '.TiV_nominal','.TiV_LDA','.TiV_QDA',
            '.ZrTi_nominal','.ZrTi_LDA','.ZrTi_QDA',
            '.TiZrYSr')

message('tectotree_all')
library(rpart)
treedata_all <- training[c(1,5:55)]
my.control <- rpart.control(xval=10, cp=0, minsplit=1)
unpruned <- rpart(AFFINITY ~ ., data=treedata_all,
                  method="class", control=my.control)
.tectotree_all <- prune(unpruned, cp=0.008)
tosave <- c(tosave,'.tectotree_all')

message('tectotree_HFS')
treedata_HFS <- get_training_data(c("AFFINITY","TiO2","La","Ce","Pr","Nd",
                                    "Sm","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                                    "Lu","Sc","Y","Zr","Nb","Hf","Ta","Pb",
                                    "Th","U","Nd143/Nd144","Sr87/Sr86",
                                    "Pb206/Pb204","Pb207/Pb204","Pb208/Pb204"))
unpruned <- rpart(AFFINITY ~ ., data=treedata_HFS,
                  method="class", control=my.control)
.tectotree_HFS <- prune(unpruned, cp=0.025)
tosave <- c(tosave,'.tectotree_HFS')

message('tectotree_ratios')
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

message('Bowen-Fenner boundary')

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
