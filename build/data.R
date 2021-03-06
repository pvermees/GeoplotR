setwd('~/Documents/Programming/GeoplotR/')
source('R/AFM.R')
source('R/tectodisc.R')
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
.AnAbOr <- IsoplotR:::fromJSON(file='inst/AnAbOr.json')
.YNb <- IsoplotR:::fromJSON(file='inst/Pearce_Y-Nb.json')
.YNbRb <- IsoplotR:::fromJSON(file='inst/Pearce_Y+Nb-Rb.json')
.YbTa <- IsoplotR:::fromJSON(file='inst/Pearce_Yb-Ta.json')
.YbTaRb <- IsoplotR:::fromJSON(file='inst/Pearce_Yb+Ta-Rb.json')
tosave <- c(tosave,'.TAS','.AnAbOr','.YNb','.YNbRb','.YbTa','.YbTaRb')

.TiZrY_nominal <- IsoplotR:::fromJSON(file='inst/TiZrY.json')
.TiZrY_LDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
attributes(.TiZrY_LDA$fit$terms)$.Environment <- NULL
.TiZrY_QDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
attributes(.TiZrY_QDA$fit$terms)$.Environment <- NULL
tosave <- c(tosave,'.TiZrY_nominal','.TiZrY_LDA','.TiZrY_QDA')

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

.BF <- construct_BF(cath)
tosave <- c(tosave,'.BF')

save(list=tosave,file="R/sysdata.rda",version=2)
