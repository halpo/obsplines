\name{fitLS}
\alias{fitLS}
\title{Fitting splines with penalized least squares.}
\description{
	Estimates the control vector for a spline fit by penalized least squares.
	The penalty being the penalty parameter times the functional inner product of the second derivative of the spline curve.
}
\usage{
fitLS(object, x, y, penalty = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ The \code{SplineBasis} object to be used to make the fit }
  \item{x}{ predictor variable. }
  \item{y}{ response variable. }
  \item{penalty}{ The penalty multiplier. }
}
\details{
	For numeric vector y, and x, and a set of basis functions, represented in \code{object}, defined on the knots \eqn{(k_0,\ldots,k_m)}.
	The likelihood is defined by\deqn{
	\sum\limits_{i=1}^n(y_i-b(x_i)\mu) + \int\limits_{k_0}^{k_m} \mu^Tb^{\prime\prime}(t)^Tb^{\prime\prime}(t)\mu dt
	}{sum_i (y_i-b(x_i)mu) + integral mu^T b''(t)^T b''(t) mu dt}
	The function estimates \eqn{\mu}.
}
\value{
  a vector of the control points.
}

\seealso{ \code{\linkS4class{SplineBasis}} }
\examples{
knots<-c(0,0,0,0:5,5,5,5)
base<-SplineBasis(knots)
x<-seq(0,5,by=.5)
y<-exp(x)+rnorm(length(x),sd=5)
fitLS(base,x,y)
}
\keyword{smooth}

