\name{SpineBasis}
\alias{SplineBasis}
\alias{OrthogonalizedSplineBasis}
\alias{evaluate}
\alias{orthogonalize}
\alias{EvaluateBasis}
\alias{print.SplineBasis}
\alias{plot.SplineBasis}
\alias{dim.SplineBasis}
\alias{OuterProdSecondDerivative}
\alias{OrthogonalSplineBasis}
\alias{OBasis}
\alias{OrthogonalizeBasis}
  
\title{Creating and Using SplineBasis Objects.}
\description{
	Show  how to from the \code{SplineBasis} objects and effectivly use them.
}
\usage{
	SplineBasis(knots, order=4)
	OrthogonalSplineBasis(knots, order=4)
	OBasis(...)
	orthogonalize(object,...)
	OrthogonalizeBasis(object,...)
	dim.SplineBasis(x)
	plot.SplineBasis(x,y,...)
	evaluate(object, t, ...)
}
\arguments{
	\item{knots}{The knots used to define the basis functions.}
	\item{order}{Order of the spline fit.(degree= order-1)}
	\item{object}{A SplineBasis object}
	\item{x}{either a SplineBasis or an OrthogonalSplineBasis object.}
	\item{y}{Should not be used}
	\item{t}{Points for evaluation.}
	\item{...}{Other arguments either ignored or passed onto other functions.}
}
\details{
	\code{SplineBasis} produces an object representing the basis functions used in spline fitting.  Provides a compact easily evaluated representation of the functions.  Produces a class of object \code{SplineBasis}.
	
	\code{orthogonalize} and \code{OrthogonalizeBasis} will provide an exact orthogonalization of the spline basis functions in \code{object}.  orthogonalize is an S4 method wrapper for OrthogonalizeBasis. Produces an object of class \code{OrthogonalSplineBasis}.
	
	\code{OrthogonalSplineBasis} is a shortcut to obtain a set of orthogonalized basis functions from the knots. \code{OBasis} is an alias for \code{OrthogonalSplineBasisA}. Both provide an object of class \code{OrthogonalSplineBasis}.  The class \code{OrthogonalSplineBasis} inherits directly from \code{SplineBasis} meaning all functions that apply to SplineBasis functions also apply to the orthogonalized version.
	
	Using \code{dim} on a SplineBasis object gives the dimentions of the internal storrage array. Of particular use is \code{dim(x)[2]} that will give the number of functions in the set.  Results in a vector of length 3.
	
	Plotting a \code{SplineBasis} object plots the basis functions.
	
	\code{evaluate} will evaluate the basis functions at the point defined in \code{x}.  Value is a p by n array, where p is the number of functions in the basis, and n is the length of \code{t}.
}
\value{
 See Details.
}
\seealso{
\code{\link{spline}}, \link{pfda-package}{pfda}
}
\references{
\emph{General matrix representations for Bsplines} KaihuaiQin, The Visual Computer 2000 16:177–186
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
