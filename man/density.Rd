\name{density_lim}
\alias{density_lim}
\alias{domain}
\alias{cutoff}
\alias{fit.density}
\title{ Analyze the eigenvalue spectrum of random matrices }
\description{
Supporting functions for analyzing the eigenvalue spectrum of random matrices
}
\usage{
density_lim(...)
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
  model <- create(WishartModel,20,50, sigma=1)
  plot(density_lim(model))

  domain(model)

  cutoff(rmatrix(model))
}
\keyword{ array }
