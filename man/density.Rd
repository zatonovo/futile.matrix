\name{dmatrix}
\alias{dmatrix}
\alias{domain}
\alias{cutoff}
\alias{fit.density}
\title{ Analyze the eigenvalue spectrum of random matrices }
\description{
Supporting functions for analyzing the eigenvalue spectrum of random matrices
}
\usage{
dmatrix(...)
domain(...)
cutoff(...)
fit.density(...)
}
\arguments{
  \item{\dots}{ Function arguments. See details for info. }
}
\details{
Calculate the upper bound of the noise spectrum.
}
\value{
'cutoff' returns the eigenvalue associated with the upper bound of the noise
spectrum.
}
\author{ Brian Lee Yung Rowe }
\examples{
  model <- create(WishartModel,20,50, sd=1)
  bounds <- domain(model)
  x <- seq(bounds[1], bounds[2], length.out=50)
  plot(x,dmatrix(x, model))

  cutoff(rmatrix(model))
}
\keyword{ array }
