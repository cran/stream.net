\name{net.object}
\alias{net.object}
\title{Data Structure for Stream Network Object}
\description{
 Components of the \code{\link{net.object}}.
}
\format{
 The network data structure, as returned from \code{net.qmodel}
 and \code{net.arcinput}, is a list with three components:

 \tabular{ll}{
   links    \tab a data frame of the topological structure\cr
   segs     \tab a data frame of the segment structure (pieces of links)\cr
   cords    \tab a list of the x and y coordinates for each segment
 }

 links has the components
   \tabular{ll}{
     lid      \tab identifier for the link\cr
     parent   \tab index of parent link (or 0 if root)\cr
     left     \tab index of left child (or 0 if terminal)\cr
     right    \tab index of right child (or 0 if terminal)\cr
     depth    \tab topological depth of link (number of links from root)\cr
     first    \tab index of first segment of this link\cr
     last     \tab index of last segment of this link\cr
     strahler \tab Strahler order\cr
     shreve   \tab Shreve order
   }

 segs has the components
   \tabular{ll}{
     sid      \tab identifier for the segment (for external ref)\cr
     link     \tab index to parent link\cr
     nxt      \tab index of the next segment for this link\cr
     prev     \tab index of the previous segment for this link\cr
     up       \tab if segment (and coordinates) oriented upstream = 1, else = 0\cr
     length   \tab segment coordinate length
   }

 cords has the components
   \tabular{ll}{
     x        \tab x coordinates with a sublist for each segment\cr
     y        \tab y coordinates with a sublist for each segment
   }
 (The names for the \code{x} and \code{y} sublists are
 the respective \code{sid}s.)
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.arcinput}}
  \code{\link{net.qmodel}}
}
\examples{
  net <- net.qmodel (5)
  net
}
\keyword{data}
