\name{peek}
\alias{peek}
\title{ Peeks inside a matrix or vector }
\description{
  Peek is a simple utility to conveniently look at a portion of a matrix. This
  is similar to head and tail but provides a 2-dimensional slice instead of a
  complete row.
}
\usage{
peek(x, upper = 5, lower = 1)
}
\arguments{
  \item{x}{ Any object that supports subsetting }
  \item{upper}{ The upper bound in the subsetting }
  \item{lower}{ The lower bound in the subsetting }
}
\value{
  Peek returns a subset of the original matrix, data.frame, etc.
}
\author{ Brian Lee Yung Rowe }
\examples{
  m <- matrix(c(1,3,4,2, 5,10,11,2, 3,42,8,22, 23,15,3,8), ncol=4)
  peek(m, 2)
}
\keyword{ array }
