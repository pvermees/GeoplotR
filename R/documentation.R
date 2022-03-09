#' @title training data
#' @description MOIB, OIB and IAB training set from Vermeesch (2006)
#' @name training
#' @docType data
#' @keywords data
#' @references Vermeesch, P., 2006, Tectonic discrimination diagrams
#'     revisited: Geochemistry, Geophysics, and Geosystems, 7, Q06017,
#'     doi:10.1029/2005GC001092
#' @examples
#' data(training,package='GeoplotR')
#' TiZrY(wtpct2ppm(training[,'TiO2']),training[,'Zr'],Y=training[,'Y'])
NULL

#' @title test data
#' @description MOIB, OIB and IAB test set from Vermeesch (2006)
#' @name test
#' @docType data
#' @keywords data
#' @references Vermeesch, P., 2006, Tectonic discrimination diagrams
#'     revisited: Geochemistry, Geophysics, and Geosystems, 7, Q06017,
#'     doi:10.1029/2005GC001092
#' @examples
#' data(test,package='GeoplotR')
#' TiZrY(wtpct2ppm(test[,'TiO2']),test[,'Zr'],test[,'Y'])
NULL

#' @title A-F-M data
#' @description Calc-alkaline and tholeiitic data from the Cascade
#'     Mountains and Iceland
#' @name cath
#' @docType data
#' @keywords data
#' @references Vermeesch, P. and Pease, V. A genetic classification of
#'     the tholeiitic and calc-alkaline magma series, Geochemical
#'     Perspective Letters (in review).
#' @examples
#' data(cath,package='GeoplotR')
#' ca <- cath[cath$affinity=='ca',]
#' A <- cath[,'Na2O']+cath[,'K2O']
#' F <- cath[,'FeOT']
#' M <- cath[,'MgO']
#' AFM(A,F,M,pch=16,col=cath$affinity)
#' legend('topleft',legend=c('Cascades','Iceland'),pch=16,col=c(1,2))
NULL
