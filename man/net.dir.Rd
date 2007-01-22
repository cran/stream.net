\name{net.dir}
\alias{net.dir}
\title{Calculate Direction Matrix for Stream Segments}
\description{
 Compute the upstream/downstream/mixed direction matrix
 for stream segments .
}
\usage{
  net.dir (dist)
}
\arguments{
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}.}
}
\details{
 Directions are +1 for upstream, -1 for downstream, and 0 for
 both.  The matrix is skew-symmetric and is read from rows to
 columns, i.e., row 1, column 3 is +1 if segment 3 is upstream
 of segment 1 and -1 if downstream.  Upstream and downstream
 are defined as being in direct linkage; any pair that have
 mixed upstream and downstream directions have direction 0.
}
\value{A square matrix of dimension the number of segments.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.dist}}
}
\examples{
  # Q model random net
  net <- net.qmodel (10)

  # add segments
  net <- net.addsegs (net)

  # distance matrix
  dmat <- net.dist (net, ends=0.5, method="coordinate")

  # directions
  tdmat <- net.dir (dmat)
  net.map (net, segatt="sid")
  table (tdmat[5,]) # includes self in zero count
}
\keyword{manip}
