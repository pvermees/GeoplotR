% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Vermeesch2006.R
\name{TiZrY}
\alias{TiZrY}
\title{Ti-Zr-Y}
\usage{
TiZrY(
  Ti = NULL,
  Zr = NULL,
  Y = NULL,
  type = c("Pearce", "LDA", "QDA"),
  ternary = TRUE,
  pch = 21,
  bg = NULL,
  show.labels = FALSE,
  short = TRUE,
  ...
)
}
\arguments{
\item{Ti}{vector with Ti concentrations (ppm)}

\item{Zr}{vector with Zr concentrations (ppm)}

\item{Y}{vector with Y concentrations (ppm)}

\item{type}{either \code{'LDA'} for linear discriminant analysis,
\code{'QDA'} for quadratic discriminant analysis, or
\code{'Pearce'} for the nominal decision boundaries of Pearce
and Cann (1973).}

\item{ternary}{logical. If \code{FALSE}, produces a logratio plot.
Only used if \code{type} is \code{LDA} or \code{QDA}.}

\item{pch}{plot character. See \code{?par} for details.}

\item{bg}{fill colour for the plot character.}

\item{show.labels}{logical. If \code{TRUE}, labels the
discrimination fields.}

\item{short}{logical. If \code{TRUE}, uses abbreviated labels for
the discrimination fields.}

\item{...}{additional arguments for the generic \code{points}
function.}
}
\value{
if \code{type='LDA'} or \code{type='QDA'}, a list with
    components \code{class}, \code{posterior} and \code{x};
    otherwise a table with labels for \code{MORB}, \code{IAB} and
    \code{OIB}.
}
\description{
Ti-Zr-Y tectonic discrimination diagram
}
\examples{
data(test,package='GeoplotR')
TiZrY(Ti=wtpct2ppm(test[,'TiO2']),
      Zr=test[,'Zr'],Y=test[,'Y'],
      type='QDA',plot='ternary')
}
\references{
Pearce, J. A., and Cann, J. R., 1973, Tectonic setting
    of basic volcanic rocks determined using trace element
    analyses: Earth and Planetary Science Letters, v. 19, no. 2,
    p. 290-300.
}
