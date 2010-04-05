\name{read.matrix}
\alias{read.matrix}
\alias{assignMatrixSparse}
\alias{assignMatrixTriplet}
\alias{assignMatrixDense}
\title{ Read a sparse matrix from a file and return a matrix }
\description{
  Reading matrices from files can be time consuming depending on the size of
  the matrix. 'read.matrix' implements a fairly efficient routine to read in 
  sparse matrices and returning dense matrix counterparts.
}
\usage{
read.matrix(file, header = FALSE, skip = 1, row.ids = NULL, col.ids = NULL, 
  colClasses = c('character','character','numeric'),
  assign.fn = assignMatrixTriplet, filter.fn = NULL, ...)
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
  \item{assign.fn}{ The function to use to construct the sparse representation
    that is then converted to a dense matrix }
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
  leverages optimized conversions in the slam and SparseM packages.

  The primary benefit of this function is that named rows and columns can be
  used as opposed to integer indexes, as compared to the \code{slam} package.
  The other main motivation is that if the memory is available, dense matrix
  calculations can be faster than their sparse counterparts, not to mention
  having a wider range of operators available.

  When header == TRUE, the row names and/or column names are read from the
  file. The names are expected to be comma separated in a single line.

  Various methods can be used to construct a sparse matrix representation
  that is used as the basis for constructing the dense matrix. The default
  function, 'assignMatrixTriplet' uses the 'slam' package under the hood. The
  'assignMatrixSparse' uses SparseM and a CSR representation, while 
  'assignMatrixDense' is a naive brute force construction that is not optimized.
}
\value{
  A matrix object generated from sparse data.
}
\author{ Brian Lee Yung Rowe }
\examples{
}
\keyword{ array }
