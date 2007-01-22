\name{net.maxupslope}
\alias{net.maxupslope}
\title{Calculate Maximum Upslope Matrix for Stream Segments}
\description{
 Compute the asymmetric maximum upslope matrix for stream
 segments.
}
\usage{
  net.maxupslope (net, slopes)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{slopes}{vector of "slopes" in segment order.}
}
\details{
 Uses the distance matrix algorithm to calculate the maximum 
 upstream slope between two segments.  The maximum is the 
 maximum of the slopes of all segments between the two 
 segments, not including either of them.  The matrix is read
 from row to column, i.e., the maximum slope from segment 1 to
 segment 3 is in row 1, column 3.

 The function was written for application to topographic
 slope but any other attribute could be used of course.
}
\value{A square matrix of dimension the number of segments.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.dist}}
}
\examples{
  # Q model random net
  net <- net.qmodel (5)

  # add segments
  net <- net.addsegs (net)

  # upslope distance matrix
  slopes <- max (net$links$shreve) - 
    net$links$shreve[net$segs$link]
  slopes <- sapply (slopes + round (runif (nrow (net$segs)), 2),
    function (x) max (0, x))
  net.map (net, segatt=slopes)
  slopemat <- net.maxupslope (net, slopes=slopes)
}
\keyword{manip}
