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

}
\value{
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- create(WignerModel, 10)
  m <- rmatrix(model)
}
\keyword{ array }
