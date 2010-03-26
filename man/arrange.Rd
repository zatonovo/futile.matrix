\name{arrange}
\alias{arrange}
\title{ Order matrix elements based on rownames and/or colnames }
\description{
  'arrange' is a convenience function for ordering rows and columns within a
  matrix.
}
\usage{
arrange(m, order.rows = TRUE, order.cols = TRUE, comparator = NULL)
}
\arguments{
  \item{m}{ A matrix whose rows and/or columns are to be arranged }
  \item{order.rows}{ Whether rows are to be ordered. Defaults to TRUE }
  \item{order.cols}{ Whether columns are to be ordered. Defaults to TRUE }
  \item{comparator}{ A function to define the ordering. Currently unused }
}
\details{
  To ensure proper operations are performed, the ordering of rows and columns
  within matrix operands must be consistent. Arrange conveniently performs
  this ordering.

  In the future, a comparator will be added so that custom orderings can be
  applied to the function.
}
\value{
  A matrix whose rows and/or columns have been ordered. By default, both rows
  and columns are ordered.
}
\author{ Brian Lee Yung Rowe }
\examples{
  library(datasets)
  arrange(state.x77)
}
\keyword{ array }
