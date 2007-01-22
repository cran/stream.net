\name{marys.aat}
\docType{data}
\alias{marys.aat}
\title{Topology Data for Marys River, Oregon}
\description{
 The "from" and "to" pointers for each arc in the Marys River
 Arc/Info coverage.  The Marys River drains an 800 square kilometer
 watershed in western Oregon.
}
\usage{marys.aat}
\format{
 This is an extraction from the Arc Attribute Table (AAT) of the
 coverage of the Marys River derived from the source below.  The
 format is comma separated values (.csv) and the three columns
 are, as stipulated in \code{\link{net.arcinput}}:

 \tabular{ll}{
   arcid  \tab integer Arc/Info User-Id for each arc (segment)\cr
   from   \tab pointer to the preceding arc (either up or downstream)\cr
   to     \tab pointer to the succeeding arc (likewise)
 }
}
\source{US Geological Survey 1:100,000 scale Digital Line Graph}
\seealso{
  \code{\link{net.arcinput}}
  \code{\link{marys.lin}}
}
\keyword{datasets}
