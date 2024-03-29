% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{xyztest}
\alias{xyztest}
\title{xyz testing function}
\usage{
xyztest(
  fname,
  xyzlab = NULL,
  show.labels = TRUE,
  short = FALSE,
  polygons = FALSE,
  smooth = FALSE,
  ...
)
}
\arguments{
\item{fname}{path to a json file}

\item{xyzlab}{(optional) 3-element vector of corner labels.}

\item{show.labels}{logical. If \code{TRUE}, labels the decision
fields.}

\item{short}{logical. If \code{TRUE}, uses short labels.}

\item{polygons}{logical. If \code{TRUE} plots the decision
polygons. If \code{TRUE}, plots the decision lines}

\item{smooth}{logical. If \code{TRUE}, plots lines as
\code{xspline}s.}

\item{...}{additional arguments to the generic plot function}
}
\value{
a list of text strings
}
\description{
tests json files with ternary decision boundaries
}
\examples{
fn <- system.file('TiZrY.json',package='GeoplotR')
xyztest(fn)
}
