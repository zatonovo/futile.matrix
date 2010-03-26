\name{expand}
\alias{expand}
\title{ Expand a matrix to larger dimensions, filling in new entries }
\description{
  Provides a convenient mechanism for expanding the dimensions of a matrix with
  a specified default value. This is particularly useful if the matrix needs to
  match dimensions with another matrix. Expand can take either another matrix
  or a set of dimnames to grow the matrix. In either case, the expanded matrix
  will have dimensions that match the explicitly or implicitly specified 
  dimnames.
}
\usage{
expand(m, target, default = 0)
}
\arguments{
  \item{m}{ A matrix to expand }
  \item{target}{ A list containing the dimnames or a matrix that contains dimnames }
  \item{default}{ The default value to use for the new entries }
}
\details{
  To properly expand m to target, the rownames and colnames of m must be a
  strict subset of target's rownames and colnames. If this requirement is not
  satisfied, the behavior is currently undefined (although most likely an error
  will result). In the future, the behavior might be configurable to drop those
  rows/columns that are not in target's rownames/colnames.

  In general, 'expand' tries to err on the side of accomodation, although the
  implementation is incomplete. If target is a list, then the format is the 
  same as when constructing a matrix and passing dimnames as an argument. 
  Currently, only a list or a matrix are supported. This could be relaxed in 
  the future to any object that has a rownames and colnames.

  Note that a current limitation/feature in expand is that it orders the 
  resulting matrix by rows and columns. More precise control needs to be 
  provided here, with the default being the ordering of the rows and columns
  conforming to target.

}
\value{
  The expanded matrix with rows and columns corresponding to the target 
  dimnames.}
\author{ Brian Lee Yung Rowe }
\examples{
  rows.m <- c('row.1', 'row.2', 'row.3')
  cols.m <- c('col.1', 'col.2')
  rows.n <- c(rows.m, 'row.4')
  cols.n <- c(cols.m, 'col.3')
  m <- matrix(c(1,4,7,2,5,8), ncol=2, dimnames=list(rows.m,cols.m))
  n <- matrix(c(1,4,7,10,2,5,8,11,3,6,9,12), ncol=3, 
    dimnames=list(rows.n,cols.n))
  expand(m, n)
}
\keyword{ array }
