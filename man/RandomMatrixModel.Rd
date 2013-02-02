\name{RandomMatrixModel}
\alias{RandomMatrixModel}
\alias{JacobiModel}
\alias{WignerModel}
\alias{WishartModel}
\alias{Ensemble}
\alias{print.Ensemble}
\title{ Type constructors for random matrices and ensembles of random matrices }
\description{
  Provides type constructors for creating random matrices. Various
  studies can be initiated afterward.
}
\usage{
RandomMatrixModel(...)
WignerModel(...)
WishartModel(...)
JacobiModel(...)
Ensemble(...)
\method{print}{Ensemble}(x, ...)
}
\arguments{
  \item{x}{An Ensemble to print}
  \item{\dots}{Additional arguments}
}
\details{
RandomMatrixModel(real=TRUE, ...)

WignerModel(n, ...)

WishartModel(n, m, sd=1, ...)

JacobiModel(n, m1, m2, ...)

Ensemble(count, model)

real - Whether the matrix has real components or not

n - Number of rows

m - Number of columns

m1 - Number of columns

m2 - Number of columns

sd - Standard deviation of the sample population

count - Number of matrices in the ensemble

model - The random matrix model to use
}
\value{
Returns a model type. Use with 'rmatrix' or 'ensemble' to generate actual 
matrices.
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- WignerModel(10)
  m <- rmatrix(model)

  e <- Ensemble(20, model)
}
\keyword{ array }
