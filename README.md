# GeoplotR

**GeoplotR** is an **R** package to classify geochemical data.

## Web interface

**GeoplotR** can be run online, offline, and from the command
line. The easiest way to get started with the program is via its web
interface:

[http://isoplotr.es.ucl.ac.uk/geoplotr/](http://isoplotr.es.ucl.ac.uk/geoplotr/)

## Installation for offline use

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
remotes::install_github('pvermees/GeoplotRgui')
```

Note that the command line functionality and graphical user interface
(GUI) are stored in two separate R packages. To start the GUI on
your own computer:

```
GeoplotRgui::GeoplotRgui()
```

## Command line

Advanced users can get access to all the functions in the GUI, and
more, from the R console. To view the contents of the API:

```
help(package="GeoplotR")
```

Here are some command line examples:

1. Quadratic discriminant analysis of Ti-Zr-Y data:

   ```
   library(GeoplotR)
   data(test,package='GeoplotR')
   TiZrY(Ti=wtpct2ppm(test[,'TiO2'],oxide='TiO2'),
   Zr=test[,'Zr'],Y=test[,'Y'],type='QDA',ternary=TRUE)
   ```

2. Total-Alkali-Silica diagram:

   ```
   out <- TAS(Na2O=test[,'Na2O'],K2O=test[,'K2O'],SiO2=test[,'SiO2'])
   ```

   Where `out` lists the classification of each sample in `test`.  

3. Plot basalt compositions from Iceland and the Cascade Mountains on
an A-F-M diagram, and mark the two-stage calc-alkaline/tholeiitic
decision boundary of Vermeesch and Pease (Geochemical Perspective
Letters, 2021) on it:

   ```
   data(cath,package='GeoplotR')
   A <- cath[,'Na2O']+cath[,'K2O']
   F <- cath[,'FeOT']
   M <- cath[,'MgO']
   AFM(A,F,M,ternary=TRUE,decision=TRUE,
   twostage=TRUE,pch=21,bg=cath[,'affinity'])
   legend('topleft',legend=c('Cascades','Iceland'),
   pch=21,pt.bg=c(1,2))
   ```

4. Plot the Iceland data in logratio space and visualise their
*Bowen-Fenner* indices as a kernel density estimate:

   ```
   data(cath,package='GeoplotR')
   TH <- (cath[,'affinity']=='th')
   A <- cath[TH,'Na2O']+cath[TH,'K2O']
   F <- cath[TH,'FeOT']
   M <- cath[TH,'MgO']
   AFM(A,F,M,ternary=FALSE,decision=TRUE,twostage=TRUE)
   ```

## Author

[Pieter Vermeesch](http://ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License
