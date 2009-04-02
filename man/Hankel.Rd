\name{Hankel}
\alias{Hankel}
\title{ Generating a Hankel Matrix }
\description{
  Functions to generate a Hankel matrix.
}
\usage{
Hankel(x, nrow = length(x), ncol = length(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ numeric vector to specify the entries of the matrix.  Should have an even number of entries.}
  \item{nrow}{ integer, must be at most length(x) }
  \item{ncol}{ integer, must be at most length(x) }
}
\details{
	Computes a Hankel matrix.  If we denote the vector $x=(x_1,\dots,x_n)$ the Hankel matrix is defined and formed as\[
	H=\left(
\begin{array}{ccccc}
	x_1&x_2&x_3&\cdots&x_{1/2}\\
	x_2&x_3&&\adots&\vdots\\
	x_3&&\adots&&\vdots\\
	\vdots&\adots&&&\vdots\\
	x_{1/2}&\cdots&\cdots&\cdots&x_n
\end{array}\right).
	\]
}
\value{
  a matrix as defined above.
}
\references{ ~put references to the literature/web site here ~ }
\author{Andrew Redd <aredd at stat.tamu.edu>}

% \seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
Hankel(1:6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }
