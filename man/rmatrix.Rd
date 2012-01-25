\name{rmatrix}
\alias{rmatrix}
\alias{rcomp}
\alias{ct}
\alias{hermitian}
\alias{symmetric}
\alias{eigenvalues}
\alias{max_eigen}
\title{ Operate on random matrices }
\description{
  Functions for generating random matrices. Also create ensembles of random
  matrices for analysis.
}
\usage{
rmatrix(...)
rcomp(...)
ct(...)
hermitian(...)
symmetric(...)
eigenvalues(...)
max_eigen(...)
}
\arguments{
  \item{\dots}{ Function arguments. See details for info. }
}
\details{
rmatrix(model)
'rmatrix' generates a random matrix. The model argument must be of type
RandomMatrixModel. Numerous sub-types (e.g. WignerModel, WishartModel) are
supported generating the appropriate type of random matrix.

'rcomp' generates a random complex number

'eigenvalues' is a convenience function that returns matrix eigenvalues.

'max_eigen' operates on an Ensemble and returns the maximum eigenvalue for each
matrix in the ensemble.

'hermitian' and 'symmetric' are convenience functions to generate these types
of random matrices. Not fully implemented.
}
\value{
  'rmatrix' returns a random matrix based on the specified model.

  'eigenvalues' returns a vector of the eigenvalues.
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- create(WignerModel, 10)
  m <- rmatrix(model)

  e <- create(Ensemble, 20, model)
  hist(max_eigen(e), freq=FALSE)
}
\keyword{ array }
