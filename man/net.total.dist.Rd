\name{net.total.dist}
\alias{net.total.dist}
\title{Calculate Distance Matrix for Stream Segments}
\description{
 Compute the symmetric total distance matrix for stream
 segments as the sum of upstream and downstream distances.
}
\usage{
  net.total.dist (dist)
}
\arguments{
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}.}
}
\value{A square matrix of dimension the number of segments.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.dist}}
  \code{\link{net.dir}}
  \code{\link{net.prox}}
}
\examples{
  # Q model random net
  net <- net.qmodel (10)

  # add segments
  net <- net.addsegs (net)

  # distance matrix
  net.map (net, segatt=round (net$segs$length, 2))
  dmat <- net.dist (net, ends=0.5, method="coordinate")

  # total distances
  tdmat <- net.total.dist (dmat)
  hist (tdmat, col="gray", main="Total Distances")
}
\keyword{manip}
