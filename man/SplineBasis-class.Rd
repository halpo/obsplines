\name{SplineBasis-class}
\docType{class}
\alias{SplineBasis-class}
\alias{OrthogonalSplineBasis-class}
\alias{deriv,SplineBasis-method}
\alias{dim,SplineBasis-method}
\alias{evaluate,SplineBasis,numeric-method}
\alias{integrate,SplineBasis-method}
\alias{orthogonalize,SplineBasis-method}
\alias{plot,SplineBasis,missing-method}

\title{Classes "SplineBasis" and "OrthogonalSplineBasis"}
\description{Contains the matrix representation for spline basis functions.  The OrthongonalSplineBasis class has the basis functions orthogonalized.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{SplineBasis}(knots, order)} or to generate orthogonal spline basis functions directly \code{\link{OrthogonalSplineBasis}(knots, order)} or the short version \code{\link{OBasis}(knots,order)}.
}
\section{Slots}{
	 \describe{
    \item{\code{knots}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{order}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Matrices}:}{Object of class \code{"array"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{deriv}{\code{signature(expr = "SplineBasis")}: Computes the derivative of the basis functions.  Returns an object of class SplineBasis.}
    \item{dim}{\code{signature(x = "SplineBasis")}: gives the dim as the order and number of basis functions.  Returns numeric of length 2.}
    \item{evaluate}{\code{signature(object = "SplineBasis", x = "numeric")}: Evaluates the basis functions and the points provided in x.  Returns a matrix with \code{length(x)} rows and \code{dim(object)[2]} columns.}
    \item{integrate}{\code{signature(object = "SplineBasis")}: computes the integral of the basis functions defined by \[\int\limits_{k_0}^x b(t)dt\] where $k_0$ is the first knot.  Returns an object of class SplineBasis.}
    \item{orthogonalize}{\code{signature(object = "SplineBasis")}: Takes in a SplinesBasis object, computes the orthogonalization transformation and returns an object of class OrthogonalSplineBasis.}
    \item{plot}{\code{signature(x = "SplineBasis", y = "missing")}: Takes an object of class SplineBasis and plots the basis functions for the domain defined by the knots in object.}
	 }
}
\references{
\emph{General matrix representations for Bsplines} KaihuaiQin, The Visual Computer 2000 16:177–186
}
\author{Andrew Redd <aredd at stat.tamu.edu>}

\seealso{
\code{\link{SplineBasis}}
}
\examples{
showClass("SplineBasis")

knots<-c(0,0,0,0:5,5,5,5)
(base <-SplineBasis(knots))
(obase<-OBasis(knots))
plot(base)
plot(obase)
}
\keyword{classes}
\keyword{hplot}
