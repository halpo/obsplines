\name{OrthogonalizeBasis}
\alias{OrthogonalizeBasis}
\title{ Orthogonalize a Spline Basis }
\description{
  Specific function for orthogonalizing the functions in a  SplineBasis object.
}
\usage{
OrthogonalizeBasis(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A SplineBasis object}
  \item{\dots}{ ignored }
}
\value{
  An \code{\linkS4class{OrthogonalSplineBasis}} object.
}
\author{ Andrew Redd <aredd at stat.tamu.edu> }
\seealso{ \code{\linkS4class{OrthogonalSplineBasis}}, \code{\linkS4class{SplineBasis}},\code{\link{orthogonalize}} }
\keyword{ classes }
