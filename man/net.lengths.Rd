\name{net.lengths}
\alias{net.lengths}
\title{Calculate Stream Segment Lengths}
\description{
 Calculate lengths of segments for a \code{\link{net.object}}
 using distance between all coordinate points along segment.
}
\usage{
  net.lengths (net)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
}
\details{
 Intended primarily for internal use. A \code{\link{net.object}}
 created by \code{\link{net.arcinput}} or \code{\link{net.qmodel}}
 includes the segment lengths.  Length is the sum of linear
 distances between segment coordinates.
}
\value{A vector, in segment order, of lengths.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
}
\examples{
  net <- net.qmodel (10)
  net.map (net, segatt="length")
}
\keyword{manip}
