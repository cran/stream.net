\name{net.read}
\alias{net.read}
\title{Read In a Stream Network Object}
\description{
 Read in a stream network object from the tables type,
 text format files written by \code{\link{net.write}}. 
}
\usage{
  net.read (prefix="somenet")
}
\arguments{
  \item{prefix}{character string for the prefix of
  the input file name. See details.}
}
\details{
 An alternative to \code{\link{load}} or \code{\link{source}}
 that reads table format text files.

 Read in the \code{$links}, \code{$segs}, and
 \code{$cords} of a \code{\link{net.object}} from three
 separate files that were written by \code{\link{net.write}}.
 The prefix is pasted onto suffixes \code{".links.dat"},
 \code{".segs.dat"}, and \code{".cords.dat"}, respectively,
 to make full filenames.

 Stream networks created from either \code{\link{net.qmodel}}
 or \code{\link{net.arcinput}} can be retrieved much more
 quickly using this method.
}
\value{A \code{\link{net.object}}.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.write}}
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

  # retrieve
  net <- net.read (prefix="test")
}
\keyword{IO}
\keyword{manip}
