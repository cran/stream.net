\name{net.autocorr.onelag}
\alias{net.autocorr.onelag}
\title{Calculate Autocorrelation Coefficient for Stream
 Segment Attribute}
\description{
 Compute the autocorrelation coefficient of an attribute
 for a specified separation distance (lag) between
 stream segments.
}
\usage{
  net.autocorr.onelag (net, dist, segatt, lag=1, eps=1e-6)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}, probably using \code{ends=0.5},
    and possibly using \code{method="segment"}.}
  \item{segatt}{attribute associated with each segment.
    If \code{segatt} is a character vector of length one,
    then get attribute from \code{net$segs$segatt}, else
    assume \code{segatt} is a numeric vector of length
    equal to number of segments and in correct order.}
  \item{lag}{separation distance (in \code{method} units)
    at which to calculate coefficient.}
  \item{eps}{precision of calculating lag.}
}
\details{
 Uses total distances, ignoring upstream/downstream. The
 algorithm for autocorrelation coefficients, from the
 reference below, is

 \preformatted{r(x) = sum(2 * z[i] * z[j])/sum(z[i]^2 + z[j]^2)
 }

 for all segment pairs \code{i,j} at lag \code{x}.  Values
 are in [-1, 1].
}
\value{
 A list with the following components:
   \tabular{ll}{
     r    \tab autocorrelation coefficient\cr
     num  \tab number of pairs used\cr
   }
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
  \code{\link{net.autocorr.att}}
  \code{\link{net.correlogram}}
}
\examples{
  # Q model random net
  net <- net.qmodel (5)
  net <- net.addsegs (net)
  dmat <- net.dist (net, ends=0.5, method="segment")

  # make rnorm and autocorrelated attributes
  att.rnorm <- net.addatt (net, dist="Normal",
    boundscaling=TRUE, vector=TRUE)
  att.autoc <- net.autocorr.att (net, dist=dmat, vector=TRUE,
    outscaling=TRUE) * 10

  # check autocorrelation
  net.autocorr.onelag (net, dist=dmat, segatt=att.rnorm)
  net.autocorr.onelag (net, dist=dmat, segatt=att.autoc)
}
\keyword{manip}
