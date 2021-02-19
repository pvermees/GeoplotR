setwd('~/Documents/Programming/R/GeoplotR/')
source('R/AFM.R')
source('R/tectodisc.R')
source('R/toolbox.R')

training <- read.csv('inst/training.csv',header=TRUE,check.names=FALSE)
save(training,file="data/training.rda",version=2)

test <- read.csv('inst/test.csv',header=TRUE,check.names=FALSE)
save(test,file="data/test.rda",version=2)

cath <- read.csv('inst/cath.csv',header=TRUE)
save(cath,file="data/cath.rda",version=2)

.TAS <- IsoplotR:::fromJSON(file='inst/TAS.json')
.AnAbOr <- IsoplotR:::fromJSON(file='inst/AnAbOr.json')
.AlFeTiMg <- IsoplotR:::fromJSON(file='inst/AlFeTiMg.json')
.PHT84 <- IsoplotR:::fromJSON(file='inst/PHT84.json')

am <- read.csv('inst/atomicmass.csv',header=FALSE)
.atomicmass <- am[,2]
names(.atomicmass) <- am[,1]

.oxides <- read.csv('inst/oxides.csv',header=FALSE,row.names=1)
colnames(.oxides) <- c('cation','ncat','nO')

.TiZrY_LDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
attributes(.TiZrY_LDA$fit$terms)$.Environment <- NULL

.TiZrY_QDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
attributes(.TiZrY_QDA$fit$terms)$.Environment <- NULL

library(rpart)
treedata_all <- training[c(1,5:55)]
my.control <- rpart.control(xval=10, cp=0, minsplit=1)
unpruned <- rpart(AFFINITY ~ ., data=treedata_all,
                  method="class", control=my.control)
.tectotree_all <- prune(unpruned, cp=0.008)

treedata_HFS <- get_training_data(c("AFFINITY","TiO2","La","Ce","Pr","Nd",
                                    "Sm","Gd","Tb","Dy","Ho","Er","Tm","Yb",
                                    "Lu","Sc","Y","Zr","Nb","Hf","Ta","Pb",
                                    "Th","U","Nd143/Nd144","Sr87/Sr86",
                                    "Pb206/Pb204","Pb207/Pb204","Pb208/Pb204"))
unpruned <- rpart(AFFINITY ~ ., data=treedata_HFS,
                  method="class", control=my.control)
.tectotree_HFS <- prune(unpruned, cp=0.025)

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

.BF <- construct_BF(cath)

save(.TiZrY_LDA,.TiZrY_QDA,.atomicmass,.oxides,.TAS,.AnAbOr,
     .AlFeTiMg,.PHT84,.tectotree_all,.tectotree_HFS,
     .tectotree_ratios,.BF,file="R/sysdata.rda",version=2)
