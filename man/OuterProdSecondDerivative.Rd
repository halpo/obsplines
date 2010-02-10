\name{OuterProdSecondDerivative}
\alias{OuterProdSecondDerivative}
\title{Outer Product of Second Derivatives of Spline Bases}
\description{
  Provides the functional outer product of second derivatives of a set of basis functions in a SplineBasis object.
  It a convenient form for forming a penalty on curve smoothness when fitting a spline curve.
}
\usage{
OuterProdSecondDerivative(basis)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{basis}{A SplineBasis object}
}
\value{
  A square matrix of order nrow(basis).
}
\author{ Andrew Redd <aredd at stat.tamu.edu> }
\seealso{ \code{\linkS4class{SplineBasis}},\code{\link{fitLS}} }
\keyword{ array }
\keyword{ math }
