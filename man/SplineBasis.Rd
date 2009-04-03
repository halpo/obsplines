\name{SplineBasis}
\alias{SplineBasis}
\alias{OrthogonalSplineBasis}
\alias{OBasis}

\title{Creating and SplineBasis Objects.}
\description{
	The function to create \code{SplineBasis} and \code{OrthogonalSplineBasis} Objects
}
\usage{
	SplineBasis(knots, order=4, keep.duplicates=FALSE)
	OrthogonalSplineBasis(knots, \dots)
	OBasis(\dots)
}
\arguments{
	\item{knots}{The full set of knots used to define the basis functions.}
	\item{order}{Order of the spline fit.(degree= order-1)}
	\item{keep.duplicates}{Should duplicate interior knots that could cause computation problem be kept or removed.  Defaults to false, which removes duplicate knots with a warning if duplicate interior knots are found.}
	\item{\dots}{Other arguments either ignored or passed onto other functions.}
}
\details{
	\code{SplineBasis} produces an object representing the basis functions used in spline fitting.  Provides a compact easily evaluated representation of the functions.  Produces a class of object \code{SplineBasis}.
	\code{OrthogonalSplineBasis} is a shortcut to obtain a set of orthogonalized basis functions from the knots. \code{OBasis} is an alias for \code{OrthogonalSplineBasis}. Both provide an object of class \code{OrthogonalSplineBasis}.  The class \code{OrthogonalSplineBasis} inherits directly from \code{SplineBasis} meaning all functions that apply to SplineBasis functions also apply to the orthogonalized version.
}
\value{
 Object of class SplineBasis or OrthogonalSplineBasis
}
\seealso{
\code{\linkS4class{SplineBasis}}, \code{\link{spline}}, \link{orthogonalsplinebasis-package}
}
\references{
\emph{General matrix representations for Bsplines} KaihuaiQin, The Visual Computer 2000 16:177--186
}
\examples{
knots<-c(0,0,0,0:10,10,10,10)
plot(SplineBasis(knots))
obase<-OBasis(knots)
plot(obase)
dim(obase)[2] #number of functions
evaluate(obase, 1:10-.5)
}
\keyword{hplot}
\keyword{algebra}
