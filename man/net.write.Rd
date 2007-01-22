\name{net.write}
\alias{net.write}
\title{Write Out a Stream Network Object}
\description{
 Write out a stream network object to a set of tables
 type, text format files. 
}
\usage{
  net.write (net, prefix="somenet")
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{prefix}{character string for the prefix of
  the output file name. See details.}
}
\details{
 An alternative to \code{\link{save}} or \code{\link{dump}}
 that writes table format text files.

 Writes out the \code{$links}, \code{$segs}, and
 \code{$cords} of the \code{\link{net.object}} in R table
 format to three separate files.  The prefix is pasted
 onto suffixes \code{".links.dat"}, \code{".segs.dat"},
 and \code{".cords.dat"}, respectively, to make full
 filenames.

 Stream networks created from either \code{\link{net.qmodel}}
 or \code{\link{net.arcinput}} can be saved this way for
 quicker retrieval for subsequent use.
}
\value{\code{\link{invisible}}.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.read}}
  \code{\link{net.qmodel}}
  \code{\link{net.arcinput}}
}
\examples{
  # Q model random net
  net <- net.qmodel (10)

  # add segments
  net <- net.addsegs (net)

  # save
  net.write (net, prefix="test")
}
\keyword{IO}
\keyword{manip}
