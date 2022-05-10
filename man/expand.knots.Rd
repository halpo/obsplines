\name{expand.knots}
\Rdversion{1.1}
\alias{expand.knots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expands knots for appropriate number of knots in B-splines}
\description{
This function is for convenience of specifying knots for B-splines.
Since the user usually only want to specify the interval that they are interested in the end knots are usually duplicated.
This function interprets the first and last knots as the end points and duplicates them.
}
\usage{
expand.knots(interior, order = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{interior}{The knots including all interior and endpoint knots}
  \item{order}{the order of the splines that the knots are to be used with.  Defaults to 4, being cubic splines }
}
\value{A vector of knots with the order specified as an attribute }
\author{Andrew Redd }


\seealso{
\code{\link{SplineBasis}}, ~~~
}
\examples{
(knots<-expand.knots(1:10))
plot(OBasis(knots))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
