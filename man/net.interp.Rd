\name{net.interp}
\alias{net.interp}
\title{Interpolate Sparse Attribute to all Stream Segments}
\description{
 Create an attribute for stream segments that is interpolated
 from data provided for a subset of segments.
}
\usage{
  net.interp (net, dist, samples, predict=NULL, 
    maxdist=1e32, method="inverseDistance", power=2, 
    vector=TRUE, name=NULL)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}, probably using \code{ends=0.5},
    and possibly using \code{method="segment"}.}
  \item{samples}{two column matrix or data frame with first
    column/field \code{sid}s to segments of sampled data and
    second column/field the sampled data values.}
  \item{predict}{vector of segments to which to predict,
    unless \code{NULL}, in which case predict to all other
    segments.}
  \item{maxdist}{maximum distance (in \code{method} units)
    for neighborhood.}
  \item{method}{only \code{"inverseDistance"} currently
    implemented.}
  \item{power}{exponent for \code{method="inverse distance"}
    function.}
  \item{vector}{if \code{TRUE}, return a vector else a field
    in \code{$segs}.}
  \item{name}{name of the attribute to be added
    (if \code{! vector)}).}
}
\details{
 Uses total distances, ignoring upstream/downstream.
}
\value{
 If \code{vector} is TRUE, then a vector of the attribute
 values in segment order, else a \code{\link{net.object}}
 with the attribute added to \code{$segs}.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.dist}}
}
\examples{
  # Q model random net
  net <- net.qmodel (10)
  dmat <- net.dist (net, ends=0.5, method="segment")

  # interpolation
  samples <- matrix (c(2, 4, 6, 10, 5, 1), ncol=2)
  y <- net.interp (net, dmat, samples, power=0.25)
  net.map (net, segatt=round (y, 1))
}
\keyword{manip}
