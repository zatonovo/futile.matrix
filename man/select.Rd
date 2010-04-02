\name{select}
\alias{select}
\alias{select<-}
\title{ Select a subset of a matrix and make assignments as necessary }
\description{
  'select' extracts a subset of a matrix based on regex patterns on either
  the rownames, the colnames or both. Once this subset has been selected,
  assignments can be made following standard consistency rules.
}
\usage{
select(m, row.pat = NULL, col.pat = NULL, ...)
}
\arguments{
  \item{m}{ A matrix from which to select a subset }
  \item{row.pat}{ A regular expression to use for rownames }
  \item{col.pat}{ A regular expression to use for colnames }
  \item{...}{ Additional arguments to pass to grep }
}
\details{
  Oftentimes it is useful to get at a specific subset of data within a matrix.
  In large matrices, it can be cumbersome to access specific rows and/or 
  columns using traditional subsetting methods, particularly if it is a complex
  set that is to be extracted. 'select' provides regex searching on named
  matrices to access portions of a matrix that satisfy the regex. Note that 
  'select' will work for data.frames as well.

  It is possible to assign values to the selected subset as a means to modify
  the original matrix. Standard consistency rules must be satisfied for any
  assignment operations.
}
\value{
  A matrix containing all rows and columns that satisfy the patterns given. If
  no values match, then an empty matrix will be returned.
}
\author{ Brian Lee Yung Rowe }
\examples{
  library(datasets)
  select(swiss, "Rive")

  select(swiss, col.pat="E", fixed=TRUE)

  select(swiss, row.pat='^[A-T]', col.pat="^E")

  select(swiss, "Rive", "Ed") <- min(select(swiss, "^[^R]", "Ed"))
}
\keyword{ array }
