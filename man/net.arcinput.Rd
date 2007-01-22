\name{net.arcinput}
\alias{net.arcinput}
\title{Read ESRI-Format Stream Networks}
\description{
 Read in a stream network from ESRI Arc/Info Arc Attribute
 Table (AAT) data and from Arc/Info "generate" format lines.
}
\usage{
  net.arcinput (aatname, linname, lineformat="R")
}
\arguments{
  \item{aatname}{name of AAT object; if
    \code{mode(aatname)=="character"}, then character
    string of the AAT file name, else R object with
    AAT data.}
  \item{linname}{name of lines object; if
    \code{mode(linname)=="character"}, then character
    string of the lines file name, else R object with
    lines data.}
  \item{lineformat}{\code{"arcgen"} for Arc/Info generate
    format lines; \code{"R"} for S/R format lines; only
    applicable to file input.}
}
\details{
 Arc/Info input consists of two objects, either from 
 files in text format, or existing R objects.

 First, information from the coverage AAT, unloaded from
 the tables command, and containing at least the user-id
 for each arc and the "from" and "to" node numbers.  Any
 number of arc attributes can follow.  After unloading,
 the file should be edited so that the first record, the
 "header", has the names of the fields in the file, in
 order, separated by commas.  The user-id is expected to
 have the name "arcid", the "from" node "from", and the
 "to" node "to".  Otherwise the names can be any legal R
 names (see the R reference manual).

 Second, the arc coordinates from an ungenerate command.
 The ungenerate command should use the "fixed" option (and
 the "line" option, of course).  If the lines have already
 been converted to (modified) S/R format, then argument
 \code{"lineformat"} should be set to \code{"R"}.  The
 arcids in the line file have to correspond with the arcids
 in the AAT file, of course.

 Modified S/R format lines consist of two columns, the \code{"x"}
 coordinate and the \code{"y"} coordinate.  Prior to the first
 coordinate pair of each line is a row containing the arcid in
 the \code{"x"} column and NA in the \code{"y"} column.

 It is essential that each arc have a unique user-id ("arcid"),
 that there be no nodes of degree greater than 2, that there
 be no cycles in the network (loops), and that the topology be
 clean for this process to work correctly.

 Arcs (the units in both the AAT file and the lines file)
 become segments in package 'stream.net' and links are
 determined by those arc junctions (nodes) that are of degree
 greater than two.

 Up or downstream order for arcs is determined by finding
 the mouth as the only singleton "from" arc (up order), or
 only singleton "to" arc (down order).
}
\value{A \code{\link{net.object}}.}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
  \code{\link{read.arcgenlin}}
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
\keyword{IO}
\keyword{manip}
