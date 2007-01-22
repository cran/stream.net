\name{net.map.key}
\alias{net.map.key}
\title{Draw Legend for a Map of a Stream Attribute}
\description{
 Draw a legend to accompany a map of a stream attribute.
}
\usage{
  net.map.key (x, y, labels=NULL, cex=par("cex"), 
    pch=par("pch"), size=2.5*cex, col=NULL, head="", 
    sep=0.25*cex, horizontal=FALSE, new=FALSE)
}
\arguments{
  \item{x,y}{lower left coordinates of key in
    proportional units [0, 1].}
  \item{labels}{vector of labels for classes
    if \code{NULL}, use integers 1:length(col),
    except if \code{is.null (col)} then use "1".}
  \item{cex}{\code{\link{par}} parameter for size of
    text.}
  \item{pch}{\code{\link{par}} parameter for type of
    symbols.}
  \item{size}{size of key symbols in \code{cex} units.}
  \item{col}{\code{\link{par}} parameter for color codes.
    Default ramp is yellow-red-brown.}
  \item{head}{text heading for key.}
  \item{sep}{separation in \code{cex} units between
    adjacent symbols. If \code{sep=0}, assume continuous
    scale and use \code{pch=15}, and put lables at
    breaks between squares.}
  \item{new}{if \code{TRUE}, create a new plot.}
  \item{horizontal}{if \code{TRUE}, key runs horizontal.}
}
\details{Alternative to \code{\link{legend}}.}
\value{
 A vector of the colors provided or generated returned
 as \code{\link{invisible}}.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.group}}
}
\examples{
  data (marys.aat, marys.lin, marys.pol, marys.elevslope)

  net <- net.arcinput (marys.aat, marys.lin)
  net.map (net, outline=marys.pol, col="blue")
  hist (marys.elevslope$elev, col="gray")
  hist (marys.elevslope$slope, col="gray")
  classes <- net.group (net, marys.elevslope$elev, ngroups=7,
    method="equalInterval")
  net.map (net, group=classes$group, lwd=2, outline=marys.pol)
  net.map.key (0.8, 0.1, labels=round(classes$cuts, 0), 
    sep=0, head="elev", horizontal=FALSE)
  classes <- net.group (net, marys.elevslope$slope, ngroups=7,
    method="equalInterval")
  net.map (net, group=classes$group, lwd=2, outline=marys.pol)
  net.map.key (0.8, 0.1, labels=round(classes$cuts, 2), 
    sep=0, head="slope", horizontal=FALSE)
}
\keyword{hplot}
