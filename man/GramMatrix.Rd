\name{GramMatrix}
\alias{GramMatrix}
\title{Computing the Gram Matrix for a set of Spline Basis}
\description{Function for computing the Gram matrix of a spline basis.}

\usage{
		GramMatrix(object)
}

\arguments{
  \item{object}{a \link{SplineBasis} object}
}
\details{
	Compute the Gram Matrix.  If \code{object} denotes the basis functions \eqn{b(t) = \{b_1(t),\ldots,b_J(t)\}} then the Gram Matrix is,\deqn{
		G = \int b^T(t)b(t) \mathrm{d}t
	}
	}
\value{
  a matrix as defined above.
}
\keyword{ math }
