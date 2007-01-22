\name{marys.pol}
\docType{data}
\alias{marys.pol}
\title{Boundary Coordinate Data for Marys River Watershed, Oregon}
\description{
 The projected map coordinates for the Arc/Info coverage of
 the boundary of the Marys River Watershed.  The Marys River
 drains an 800 square kilometer watershed in western Oregon.
}
\usage{marys.pol}
\format{
 These data are derived from a "generate" format representation of
 the coordinates of the single arc in the coverage of the boundary
 of the watershed of the Marys River.  The format is modified S/R
 lines in a table (.tab) and the two columns are, as stipulated in
 \code{\link{net.arcinput}}:

 \tabular{ll}{
   x  \tab x-coordinate\cr
   y  \tab y-coordinate
 }

 Modified S/R format lines consist of two columns, the \code{"x"}
 coordinate and the \code{"y"} coordinate.  Prior to the first
 coordinate pair of each line is a row containing the line or
 polygon identifier in the \code{"x"} column and NA in the
 \code{"y"} column.

 For a polygon, the first and last coordinate pairs are identical.

 The map projection for these coordinates is Universal Transverse
 Mercator, Zone 10 using the Clarke 1866 Spheroid.
}
\source{Derived from 30 meter Digital Elevation Model}
\keyword{datasets}
