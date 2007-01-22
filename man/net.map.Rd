\name{net.map}
\alias{net.map}
\title{Draw Map of Stream Network Attribute}
\description{
 Draw map of stream network by symbolizing segments
 according to a classification of segment or link
 attributes obtained from \code{\link{net.group}},
 using colors or line widths, or just according to
 a set of colors or line widths directly associated
 with the segments.
}
\usage{
  net.map (net, group=NULL, linkatt=NULL, segatt=NULL, 
    col=NULL, lwd=NULL, cex=par("cex"), new=TRUE, 
    outline=NULL, uniquegroup=FALSE)
}
\arguments{
  \item{net}{A \code{\link{net.object}}.}
  \item{group}{if not \code{NULL}, then a vector of color
    codes for each segment.  \code{names(group)} must
    match \code{sid}s.}
  \item{segatt}{if not \code{NULL} then write the
    attribute value associated with each segment in text
    format at the midpoint of the segment.
    If \code{segatt} is a character vector of length one,
    then get attribute from \code{net$segs$segatt}, else
    assume \code{segatt} is a numeric vector of length
    equal to number of segments and in correct order.}
  \item{linkatt}{if not \code{NULL} then write the
    attribute value associated with each link in text
    format at the midpoint of the link.
    If \code{linkatt} is a character vector of length one,
    then get attribute from \code{net$links$linkatt}, else
    assume \code{linkatt} is a numeric vector of length
    equal to number of links and in correct order.}
  \item{col}{\code{\link{par}} parameter for color codes
    for segments.  Default ramp is yellow-red-brown.}
  \item{lwd}{\code{\link{par}} parameter for line widths
    of segments. If \code{lwd="seq"}, generate widths in
    sequence 1:(number of groups).}
  \item{cex}{\code{\link{par}} parameter for size of text.}
  \item{new}{if \code{TRUE}, then create a new plot.}
  \item{outline}{if not \code{NULL}, then shade an outline
    outline polygon in light gray first before drawing
    network; outline format is S/R polygon format.}
  \item{uniquegroup}{if \code{TRUE}, group categories
    are not members of the integers 1:(number of groups),
    assuming \code{! is.null (group)}.}
}
\details{
 The mapping is by segment.  Link attributes can be mapped
 by grouping them with \code{\link{net.group}} where they
 will be "stretched" onto the segments.

 Text format labeling of attributes for segments or links
 is separate from (and can be in addition to) color or
 line width symbolism of attributes.

 The \code{names} of the \code{group} vector must match
 the \code{sid}s (as set, for example, in
 \code{\link{net.group}}).

 The values of \code{group} are assumed to be integers.
 If these values are not dense (continuous in 1:(number
 of groups), then set \code{uniquegroups=TRUE}.
}
\value{
 If \code{! is.null (group)} then the colors provided or
 generated else nothing in \code{\link{invisible}}.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{net.group}}
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
  net.map (net) # just draw segments in black
  n <- length (table (class1$group))
  net.map (net, group=class1$group, col="black", lwd="seq")
  n <- length (table (class2$group))
  net.map (net, group=class2$group, col=gray(0:n/n),
    linkatt="strahler")
  n <- length (table (class3$group))
  net.map (net, group=class3$group, col=rainbow(n))
}
\keyword{hplot}
