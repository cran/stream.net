\name{net.qmodel}
\alias{net.qmodel}
\title{Generate Stochastic Stream Network Object}
\description{
 Generate a stochastic network according to the Q model
 (see references).
}
\usage{
  net.qmodel (size=100, Q=0.5)
}
\arguments{
  \item{size}{magnitude of the network in number of terminal
      branches.}
  \item{Q}{probability of an internal branch rather than a 
      terminal branch (0 <= \code{Q} <= 1).}
}
\details{
 Any value of \code{Q} will generate all possible topologies
 of the given size, however, the likelihood of any topology is
 controlled by the value provided.  If \code{Q} = 0.01, then
 networks with short lengths to terminal branches are more
 likely (because branching from existing terminal branches
 is much more probable than from existing internal
 branches).  If \code{Q} = 0.99, then networks with a few long
 sequences of branches are more likely.
}
\value{A \code{\link{net.object}}.}
\references{
 Costa-Cabral MC, Burges SJ.  1997.  Sensitivity of channel
 network planform laws and the question of topologic
 randomness.  \emph{Water Resources Research}
 \bold{33}(9):2179-2197.
}
\author{Denis White, \email{white.denis@epa.gov} }
\seealso{
  \code{\link{net.object}}
}
\examples{
  net <- net.qmodel (10, Q=0.99)
  net.map (net)
}
\keyword{datagen}
