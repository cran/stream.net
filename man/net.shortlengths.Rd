\name{net.shortlengths}
\alias{net.shortlengths}
\title{Calculate Stream Segment Endpoint-only Lengths}
\description{
 Calculate lengths of segments for a
 \code{\link{net.object}} using distance between endpoints.
}
\usage{
  net.shortlengths (net)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
}
\details{
 Intended primarily for internal use. The shortlength is
 the distances between segment endpoint coordinates
 only, ignoring any intermediate points.
}
\value{A vector, in segment order, of shortlengths.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.lengths}}
}
\examples{
  data (marys.aat, marys.lin)

  net <- net.arcinput (marys.aat, marys.lin)
  summary (net.lengths (net) - net.shortlengths (net))
}
\keyword{manip}
