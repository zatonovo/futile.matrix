\name{futile.matrix-package}
\alias{futile.matrix-package}
\alias{futile.matrix}
\docType{package}
\title{
  A collection of matrix manipulation functions
}
\description{
  This is a sub package of the futile collection of R libraries. It provides
  various convenience functions for examining data within matrices as well as
  some optimized functions for reading matrices in various formats.
}
\details{
\tabular{ll}{
Package: \tab futile.matrix\cr
Type: \tab Package\cr
Version: \tab 1.2.2\cr
Date: \tab 2014-05-06\cr
License: \tab LGPL-3\cr
LazyLoad: \tab yes\cr
}
  Matrices are a widely used data type in R, which provides many fine language
  features and built-in functions for working with them. As matrices grow in 
  size, they become harder to work with, not from a language perspective, but
  from an introspection perspective. Using the subsetting functions is typically
  the approach used to extract subsets of data from a matrix, but complex
  constructions leave code unreadable as well as fragile.  The futile.matrix 
  package addresses this issue by providing subsetting functions to access 
  data in a readable way that is (hopefully) consistent with R idioms. It also
  provides functions for growing matrices in defined ways.

  While R is strong in manipulating matrices, getting matrices into R can at
  times be a challenge. Large (but not huge) matrices can be slow to load, so
  futile.matrix provides some functions that speed up the reading process. The
  implementation provides routines for reading both dense matrices and sparse
  matrices (in CSR, CSC, and triplet forms). Note that this does not attempt 
  to replace packages like bigmemory, which deals with matrices orders of 
  magnitude larger.
}
\author{
Brian Lee Yung Rowe <r@zatonovo.com>
}
\keyword{ package }
\keyword{ attribute }
\keyword{ logic }
\seealso{ \code{\link{select}}, \code{\link{expand}}, \code{\link{read.matrix}} }
\examples{
library(datasets)
select(swiss, "Rive")
select(swiss, col.pat='^E')
select(swiss, "Rive", '^E') <- -1

dimnames <- list( c(rownames(swiss), 'Zermat', 'Zurich', 'Geneva'),
  c(colnames(swiss), 'Age','Hair.Color') )
my.swiss <- expand(swiss, dimnames)
}

