# GeoplotR

**GeoplotR** is an **R** package to classify geochemical data.

## Installation

You must have **R** installed on your system (see
[http://r-project.org](http://r-project.org)).  Additionally, to
install **GeoplotR** from Github, you also need the **remotes**
package.  This can be installed by typing the following code at the R
command line prompt:

```
install.packages('remotes')
```

after which **GeoplotR** can be installed as follows:

```
remotes::install_github('pvermees/GeoplotR')
```

## Examples

Quadratic discriminant analysis of Ti-Zr-Y data:

```
library(GeoplotR)
data(test,package='GeoplotR')
TiZrY(Ti=test[,'TiO2'],Zr=test[,'Zr'],Y=test[,'Y'],type='QDA',plot='ternary')
```

Total-Alkali-Silica diagram:

```
out <- TAS(Na2O=test[,'Na2O'],K2O=test[,'K2O'],SiO2=test[,'SiO2'],pch=16)
```

Where `out` lists the classification of each sample in `test`.  

Plot basalt compositions from Iceland and the Cascade Mountains on an
A-F-M diagram:

```
data(cath,package='GeoplotR')
ca <- cath[cath$affinity=='ca',]
AFM(cath$A,cath$F,cath$M,pch=16,col=cath$affinity)
legend('topleft',legend=c('Cascades','Iceland'),pch=16,col=c(1,2))
```

Calculate the *Bowen-Fenner* index of the data:

```
bfi <- BF(cath$A,cath$F,cath$M)
oldpar <- par(mfrow=c(2,1))
hist(bfi[cath$affinity=='ca'])
hist(bfi[cath$affinity=='th'])
par(oldpar)
```

Plot the Cascades data on a logratio plot and project a kernel density
estimate of the Bowen-Fenner indices on a radial scale:

```
Cascades <- (cath$affinity=='ca')
AFM(cath$A[Cascades],cath$F[Cascades],cath$M[Cascades],ternary=FALSE,radial=TRUE)
```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
