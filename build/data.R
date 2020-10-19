setwd('~/Documents/Programming/R/GeoplotR/')
source('R/tectodisc.R')

training <- read.delim('inst/training.txt',header=TRUE,
                       sep='\t',check.names=FALSE)
save(training,file="data/training.rda",version=2)

test <- read.delim('inst/test.txt',header=TRUE,
                   sep='\t',check.names=FALSE)
save(test,file="data/test.rda",version=2)

.TAS <- IsoplotR:::fromJSON(file='inst/TAS.json')
.AbOrAn <- IsoplotR:::fromJSON(file='inst/AbOrAn.json')
.AlFeTiMg <- IsoplotR:::fromJSON(file='inst/AlFeTiMg.json')
.Pearson <- IsoplotR:::fromJSON(file='inst/Pearson.json')

am <- read.csv('inst/atomicmass.csv',header=FALSE)
.atomicmass <- am[,2]
names(.atomicmass) <- am[,1]

.TiZrY_LDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=FALSE,plot=FALSE)
attributes(.TiZrY_LDA$fit$terms)$.Environment <- NULL

.TiZrY_QDA <- construct_DA(X='Ti',Y='Zr',Z='Y',quadratic=TRUE,plot=FALSE)
attributes(.TiZrY_QDA$fit$terms)$.Environment <- NULL

save(.TiZrY_LDA,.TiZrY_QDA,.atomicmass,.TAS,.AbOrAn,
     .AlFeTiMg,.Pearson,file="R/sysdata.rda",version=2)
