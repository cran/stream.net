\name{net.addsegs}
\alias{net.addsegs}
\title{Add Segments to Links of Stream Network}
\description{
 Add segments to a \code{\link{net.object}}, according
 to either a uniform or a gamma distribution.
}
\usage{
  net.addsegs (net, dist="Uniform", min=1, max=2, 
    shape=1.35, scale=2.3)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{dist}{If \code{"Uniform"}, use uniform distribution
    for number of segments to add, else if \code{"Gamma"},
    use gamma distribution.}
  \item{min}{\code{min} argument of \code{\link{runif}}.}
  \item{max}{\code{max} argument of \code{\link{runif}}.}
  \item{shape}{\code{shape} argument of \code{\link{rgamma}}.}
  \item{scale}{\code{scale} argument of \code{\link{rgamma}}.}
}
\details{
 Add segments to a (presumably random) network.  Number
 of segments to add per link is either a random draw from
 a uniform distribution between \code{min} and \code{max},
 inclusive, or a random draw from a gamma distribution with
 parameters \code{shape} and \code{scale}.  Does not handle
 segment attributes, so should be used before attributes are
 added.  The number of segments added for any particular
 link can be zero.

 Primarily designed to be used with \code{\link{net.qmodel}}.
}
\value{A \code{\link{net.object}} with segments added.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.qmodel}}
  \code{\link{net.addatt}}
  \code{\link{runif}}
  \code{\link{rgamma}}
}
\examples{
  net <- net.qmodel (10)
  net <- net.addsegs (net)
  net.map (net, segatt="sid")
}
\keyword{datagen}
\keyword{manip}
