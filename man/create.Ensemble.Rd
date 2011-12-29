\name{create.RandomMatrixModel}
\alias{create.RandomMatrixModel}
\alias{create.JacobiModel}
\alias{create.WignerModel}
\alias{create.WishartModel}
\alias{create.Ensemble}
\title{ Type constructors for random matrices and ensembles of random matrices }
\description{
  Provides type constructors for creating random matrices. Various
  studies can be initiated afterwards.
}
\usage{
create.RandomMatrixModel(T, real=TRUE)
create.WignerModel(T, n, ...)
create.WishartModel(T, n, m, ...)
create.JacobiModel(T, n, m1, m2, ...)
create.Ensemble(T, rank, count, model)
}
\arguments{
  \item{T}{ The type }
  \item{real}{ Whether the matrix has real components or not}
  \item{n}{ Number of rows}
  \item{m}{ Number of columns}
  \item{m1}{ Number of columns}
  \item{m2}{ Number of columns}
  \item{rank}{ Rank of the matrices }
  \item{count}{ Number of matrices in the ensemble }
  \item{model}{ The random matrix model to use }
}
\value{
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- create(WignerModel, 10)
  m <- rmatrix(model)
}
\keyword{ array }
