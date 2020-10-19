# GeoR

**GeoplotR** is an **R** package to classify geochemical data.

## Prerequisites

You must have **R** installed on your system (see
[http://r-project.org](http://r-project.org)).  Additionally, to
install **GeoplotR** from Github, you also need the **remotes**
package.  This can be installed by typing the following code at the R
command line prompt:

```
install.packages('remotes')
```

## Installation

```
remotes::install_github('pvermees/GeoplotR')
```

## Examples

Quadratic discriminant analysis of Ti-Zr-Y data:

```
library(GeoplotR)
TiZrY(Ti=test[,'TIO2(WT%)'],Zr=test[,'ZR(PPM)'],Y=test[,'Y(PPM)'],
      units=c('wt%','ppm','ppm'),type='QDA',plot='ternary')
```

Total-Alkali-Silica diagram:

```
out <- TAS(Na2O=test[,'NA2O(WT%)'],K2O=test[,'K2O(WT%)'],SiO2=test[,'SIO2(WT%)'])
```

Where `out` lists the classification of each sample in `test`.

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
