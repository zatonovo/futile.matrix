\name{read.matrix}
\alias{read.matrix}
\title{ Read a sparse matrix from a file and return a matrix }
\description{
  Reading matrices from files can be time consuming depending on the size of
  the matrix. 'read.matrix' implements a fairly efficient routine to read in 
  sparse matrices and returning dense matrix counterparts.
}
\usage{
read.matrix(file, header = TRUE, format = 'triplet', filter.fn = NULL, ...)
}
\arguments{
  \item{file}{ A file or connection to read from }
  \item{header}{ Whether header lines exist defining all possible rows and
    columns. If this is false, then the defined triplet elements will produce
    the complete set of rows and columns. Defaults to TRUE }
  \item{format}{ The sparse matrix format of the file. Currently only 'triplet'
    is supported. }
  \item{filter.fn}{ An optional function used to filter/clean the input data 
    and/or row/column ids. The signature of filter.fn must have arguments for 
    data, row.ids, and col.ids}
  \item{...}{ Additional arguments to pass to the construction portion of the
    implementation }
}
\details{
  Matrices that have dimensions on the order of thousands can be slow to load 
  into R. 'read.matrix' provides an efficient implementation for reading sparse
  matrices from a file or other connection. The underlying implementation
  leverages optimized conversions in the SparseM package from CSR formatted
  matrices to dense matrices. This turns out to be faster than performing
  direct assignment on a dense matrix element by element.
}
\value{
  A matrix object generated from sparse data.
}
\author{ Brian Lee Yung Rowe }
\examples{
}
\keyword{ array }
