\name{create.RandomMatrixModel}
\alias{create.RandomMatrixModel}
\alias{create.JacobiModel}
\alias{create.WignerModel}
\alias{create.WishartModel}
\alias{create.Ensemble}
\alias{print.Ensemble}
\title{ Type constructors for random matrices and ensembles of random matrices }
\description{
  Provides type constructors for creating random matrices. Various
  studies can be initiated afterward.
}
\usage{
create.RandomMatrixModel(T, real=TRUE, ...)
create.WignerModel(T, n, ...)
create.WishartModel(T, n, m, sd=1, ...)
create.JacobiModel(T, n, m1, m2, ...)
create.Ensemble(T, count, model)
\method{print}{Ensemble}(x, ...)
}
\arguments{
  \item{T}{ The type }
  \item{real}{ Whether the matrix has real components or not}
  \item{n}{ Number of rows}
  \item{m}{ Number of columns}
  \item{m1}{ Number of columns}
  \item{m2}{ Number of columns}
  \item{sd}{ Standard deviation of the sample population}
  \item{count}{ Number of matrices in the ensemble }
  \item{model}{ The random matrix model to use }
  \item{x}{An Ensemble to print}
  \item{\dots}{ Additional arguments}
}
\value{
Returns a model type. Use with 'rmatrix' or 'ensemble' to generate actual 
matrices.
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- create(WignerModel, 10)
  m <- rmatrix(model)

  e <- create(Ensemble, 20, model)
}
\keyword{ array }
