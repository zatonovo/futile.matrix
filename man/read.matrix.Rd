\name{read.matrix}
\alias{read.matrix}
\alias{assignMatrixSparse}
\alias{assignMatrixDense}
\title{ Read a sparse matrix from a file and return a matrix }
\description{
  Reading matrices from files can be time consuming depending on the size of
  the matrix. 'read.matrix' implements a fairly efficient routine to read in 
  sparse matrices and returning dense matrix counterparts.
}
\usage{
read.matrix(file, header = FALSE, skip = 1, row.ids = NULL, col.ids = NULL, 
  colClasses = c('character','character','numeric'), filter.fn = NULL, ...)
}
\arguments{
  \item{file}{ A file or connection to read from }
  \item{header}{ Whether header lines exist defining all possible rows and
    columns. If this is false, then the defined triplet elements will produce
    the complete set of rows and columns. }
  \item{skip}{ The number of rows to skip. This assumes there is a single 
    header line, which is skipped. }
  \item{row.ids}{ If header is TRUE, the row number that defines the row.ids.
    If header == FALSE, the row.ids to use for the matrix }
  \item{col.ids}{ If header is TRUE, the col number that defines the col.ids.
    If header == FALSE, the col.ids to use for the matrix }
  \item{colClasses}{ The classes to use for the columns in the triplet file }
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

  The primary benefit of this function is that named rows and columns can be
  used as opposed to integer indexes, as compared to the \code{slam} package.

  When header == TRUE, the row names and/or column names are read from the
  file. The names are expected to be comma separated in a single line.
}
\value{
  A matrix object generated from sparse data.
}
\author{ Brian Lee Yung Rowe }
\examples{
}
\keyword{ array }
