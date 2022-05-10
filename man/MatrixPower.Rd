\name{MatrixPower}
\alias{MatrixPower}
\title{Matrix Power}
\description{
	Performs the matrix power operation.
}
\usage{
MatrixPower(A, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{ A square matrix. }
  \item{n}{ An integer telling the exponent.}
}
\details{
 Only well defined for integers the matrix power operation is a convenience function to multiply a matrix, \code{A}, with itself \code{n} times.
}
\value{
	A matrix of the same dimension as \code{A}.
}
\examples{
A<-rbind(0,cbind(diag(1:5),0))  #a nilpotent matrix
A
MatrixPower(A,3)
MatrixPower(A,5)
MatrixPower(A,6)    #Gets to a zero matrix
}
\keyword{array}
\keyword{algebra}
