\name{net.orders}
\alias{net.orders}
\title{Calculate Strahler and Shreve Stream Orders}
\description{
 Calculate Strahler and Shreve orders for a
 \code{\link{net.object}} and add to link table.
}
\usage{
  net.orders (links)
}
\arguments{
  \item{links}{link table, see \code{\link{net.object}}.}
}
\details{
 Intended primarily for internal use. A \code{\link{net.object}}
 created by \code{\link{net.arcinput}} or \code{\link{net.qmodel}}
 includes the stream orders.
}
\value{A link table.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
}
\examples{
  net <- net.qmodel (10)
  net.map (net, linkatt="strahler")
  net.map (net, linkatt="shreve")
}
\keyword{manip}
