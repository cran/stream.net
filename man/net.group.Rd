\name{net.group}
\alias{net.group}
\title{Classify Stream Attribute for Mapping}
\description{
 Classify stream segment or stream link attributes for
 use by \code{\link{net.map}}.
}
\usage{
  net.group (net, segatt=NULL, linkatt=NULL, 
    ngroups=5, method="quantile", spread=NULL)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{segatt}{attribute associated with each segment.
    If \code{segatt} is a character vector of length one,
    then get attribute from \code{net$segs$segatt}, else
    assume \code{segatt} is a numeric vector of length
    equal to number of segments and in correct order.}
  \item{linkatt}{attribute associated with each link.
    If \code{linkatt} is a character vector of length one,
    then get attribute from \code{net$links$linkatt}, else
    assume \code{linkatt} is a numeric vector of length
    equal to number of links and in correct order.}
  \item{ngroups}{the number of groups to create.}
  \item{method}{currently must be \code{"quantile"} or
    \code{"equalInterval"}.  \code{method = "quantile"}
    means quartiles, quintiles, etc., depending on ngroups.  
    \code{"method = equalInterval"} means dividing the range
    of the attribute into \code{ngroups} equal intervals.}
  \item{spread}{if not \code{NULL}, then, for
    \code{"equalInterval"} classification, set the range
    of intervals to be the two element spread vector.}
}
\details{
 One of \code{segatt} or \code{linkatt} must not be
 \code{NULL}.
}
\value{
 A list with the following components:
   \tabular{ll}{
     group  \tab vector of group numbers for each
       segment\cr
     cuts   \tab vector of break points for the
       classification
   }
 The \code{group} vector is always of length the number
 of segments.  If the attribute is provided by link, the
 attribute is associated with all segments in the link.

 The \code{names} of the \code{group} vector are set
 to the \code{sid}s so that \code{\link{net.map}} will
 work correctly.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.map}}
}
\examples{
  # Q model random net
  net <- net.qmodel (10)
  net <- net.addsegs (net)

  # classify some attributes
  class1 <- net.group (net, linkatt="shreve", 
    method="equalInterval")
  class2 <- net.group (net, segatt="length")
  att <- seq (nrow (net$segs))
  class3 <- net.group (net, segatt=att)

  # map
  n <- length (table (class1$group))
  net.map (net, group=class1$group, col=rev(terrain.colors(n)))
  n <- length (table (class2$group))
  net.map (net, group=class2$group, col=gray(0:n/n))
  n <- length (table (class3$group))
  net.map (net, group=class3$group, col=rainbow(n))
}
\keyword{manip}
