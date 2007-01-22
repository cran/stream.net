\name{net.dist}
\alias{net.dist}
\title{Calculate Distance Matrix for Stream Segments}
\description{
 Compute the asymmetric distance matrix for stream segments.
 Upstream and downstream distances are included separately.
}
\usage{
  net.dist (net, ends=0.5, method="coordinate", digits=10)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{ends}{how to handle from and to segments
    if \code{ends=0}, do not use from and to lengths
    if \code{ends=1}, use total of from and to lengths
    if \code{ends=0.5}, use half of from and to lengths.}
  \item{method}{\code{method="coordinate"} means network
    coordinate distance; \code{method="segment"} means
    distance in number of segments.}
  \item{digits}{if \code{digits=NULL}, do not round output
    matrix, else precision of rounding.}
}
\details{
 Upstream distances are found by indexing the source segment
 by its row and the destination segment by its column.
 Downstream distances are the opposite; the source is the
 column and the destination is the row.

 Algorithm adapted from that of SG Leibowitz.  This
 algorithm does produce roundoff residue and thus the
 matrix can be purged of very small non-zero values with
 rounding.
}
\value{A square matrix of dimension the number of segments.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.total.dist}}
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
  hist (dmat, col="gray", main="One Way Distances")
}
\keyword{manip}
