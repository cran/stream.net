\name{marys.lin}
\docType{data}
\alias{marys.lin}
\title{Line Coordinate Data for Marys River, Oregon}
\description{
 The projected map coordinates for the Arc/Info coverage of
 the Marys River.  The Marys River drains an 800 square kilometer
 watershed in western Oregon.
}
\usage{marys.lin}
\format{
 These data are derived from a "generate" format representation of
 the coordinates of the arcs in the coverage of the Marys River,
 that are in turn derived from the source below.  The format is
 modified S/R lines in a table (.tab) and the two columns are,
 as stipulated in \code{\link{net.arcinput}}:

 \tabular{ll}{
   x  \tab x-coordinate\cr
   y  \tab y-coordinate
 }

 Modified S/R format lines consist of two columns, the \code{"x"}
 coordinate and the \code{"y"} coordinate.  Prior to the first
 coordinate pair of each line is a row containing the arcid in
 the \code{"x"} column and NA in the \code{"y"} column.

 The map projection for these coordinates is Universal Transverse
 Mercator, Zone 10 using the Clarke 1866 Spheroid.
}
\source{US Geological Survey 1:100,000 scale Digital Line Graph}
\seealso{
  \code{\link{read.arcgenlin}}
  \code{\link{net.arcinput}}
  \code{\link{marys.aat}}
}
\keyword{datasets}
