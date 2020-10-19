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
#' TiZrY(Ti=training[,'TIO2(WT%)'],Zr=training[,'ZR(PPM)'],
#'       Y=training[,'Y(PPM)'],units=c(2,1,1))
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
#' TiZrY(Ti=test[,'TIO2(WT%)'],Zr=test[,'ZR(PPM)'],
#'       Y=test[,'Y(PPM)'],units=c(2,1,1))
NULL
