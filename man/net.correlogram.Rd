\name{net.correlogram}
\alias{net.correlogram}
\title{Calculate Correlogram for Stream Segment Attribute}
\description{
 Compute the autocorrelation coefficients of an attribute
 for a sequence of separations in distance (lags) between
 stream segments.
}
\usage{
  net.correlogram (net, dist, segatt, nlags=10)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{dist}{an upstream/downstream distance matrix from
    \code{\link{net.dist}}, probably using \code{ends=0.5}.}
  \item{segatt}{attribute associated with each segment.
    If \code{segatt} is a character vector of length one,
    then get attribute from \code{net$segs$segatt}, else
    assume \code{segatt} is a numeric vector of length
    equal to number of segments and in correct order.}
  \item{nlags}{number of lags in the output function.}
}
\details{
 Uses total distances, ignoring upstream/downstream. The
 algorithm for autocorrelation coefficients, from the
 reference below, is

 \preformatted{r(x) = sum(2 * z[i] * z[j])/sum(z[i]^2 + z[j]^2)
 }

 for all segment pairs \code{i,j} at lag \code{x}.  Values
 are in [-1, 1].

 Lag distances are in \code{method} units used in
 \code{net.dist}, and are in range \code{[0, max(dist)]}.
}
\value{
 A list with the following elements:
   \tabular{ll}{
     r    \tab vector of autocorrelation coefficients
                at each lag\cr
     lag  \tab vector of lag distances\cr
     num  \tab vector of number of pairs at each lag\cr
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
  \code{\link{net.autocorr.onelag}}
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

  # plot correlograms
  cg1 <- net.correlogram (net, dist=dmat, segatt=att.rnorm)
  cg2 <- net.correlogram (net, dist=dmat, segatt=att.autoc)
  plot (cg1$lag, cg1$r, pch=19, type="b", ylim=c(-1,1), 
    main="Random Normal")
  plot (cg2$lag, cg2$r, pch=19, type="b", ylim=c(-1,1),
    main="Autocorrelated")
}
\keyword{manip}
