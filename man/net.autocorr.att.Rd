\name{net.autocorr.att}
\alias{net.autocorr.att}
\title{Create Autocorrelated Attribute for Stream Segments}
\description{
 Create a network autocorrelated attribute for stream
 segments.
}
\usage{
  net.autocorr.att (net, dist, target=0.5, lag=1,
    outscaling=TRUE, min=0, max=1, eps=1e-6,
    vector=TRUE, name=NULL)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}, probably using \code{ends=0.5},
    and possibly using \code{method="segment"}.}
  \item{target}{goal for autocorrelation.}
  \item{lag}{separation distance (in \code{method} units)
    at which to produce autocorrelation effect.}
  \item{outscaling}{if \code{TRUE}, scale return values to
    \code{[min, max]}.}
  \item{min}{minimum value for scaling output values.}
  \item{max}{maximum value for scaling output values.}
  \item{eps}{precision of calculating lag distance.}
  \item{vector}{if \code{TRUE}, return a vector else a field
    in \code{$segs}.}
  \item{name}{name of the attribute to be added
    (if \code{! vector)}).}
}
\details{
 Uses total distances, ignoring upstream/downstream. The
 algorithm for calculating autocorrelation, from the
 reference below, is

 \preformatted{r(x) = sum(2 * z[i] * z[j])/sum(z[i]^2 + z[j]^2)
 }

 for all segment pairs \code{i,j} at lag \code{x}.  Values
 are in [-1, 1].
}
\value{
 If \code{vector} is TRUE, then a vector of the attribute
 values in segment order, else a \code{\link{net.object}}
 with the attribute added to \code{$segs}.}
}
\references{
 Henley S.  1975.  Autocorrelation coefficients from
 irregularly spaced areal data.  \emph{Computers and
 Geosciences} \bold{2}(4):437-438.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.dist}}
}
\examples{
  # Q model random net
  net <- net.qmodel (5)
  net <- net.addsegs (net)
  dmat <- net.dist (net, ends=0.5, method="segment")

  # make autocorrelated attributes
  att.autoc <- net.autocorr.att (net, dist=dmat, vector=TRUE,
    outscaling=TRUE) * 10

  # check autocorrelation
  net.autocorr.onelag (net, dist=dmat, segatt=att.autoc)
}
\keyword{manip}
