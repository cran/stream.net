\name{read.arcgenlin}
\alias{read.arcgenlin}
\title{Read In an Arc/Info Generate Format Lines File}
\description{
 Read in an Arc/Info generate format lines file and convert
 to modified S/R lines format. 
}
\usage{
  read.arcgenlin (filename, coord=c("x","y")) 
}
\arguments{
  \item{filename}{character string for the input file name.}
  \item{coord}{the names to be given to the columns for the 
    coordinates.  For stream.net, use the default.}
}
\details{
 Any Arc/Info generate format line data can be converted for
 plotting or other uses in R.  In particular,
 \code{\link{net.arcinput}} uses Arc/Info generate format lines
 to represent the coordinate structure of a stream network.

 Modified S/R format lines consist of two columns, the \code{"x"}
 coordinate and the \code{"y"} coordinate.  Prior to the first
 coordinate pair of each line is a row containing the arcid in
 the \code{"x"} column and NA in the \code{"y"} column.
}
\value{A \code{\link{net.object}}.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.arcinput}}
}
\examples{
 # No examples provided as there appears to be no way to store
 # the Arc/Info format data directly in R.
}
\keyword{IO}
\keyword{manip}
