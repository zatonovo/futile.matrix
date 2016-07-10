# Get the cutoff for the noise portion of the spectrum, which is essentially
# lambda+. Depending on the technique, there are different ways of achieving
# this.
RandomMatrixFilter(...) %as% list(...)
MaximumLikelihoodFit(...) %as% list(...)

#' Calculate the upper bound of the noise spectrum
#'
#' Finds the cutoff point between the noise portion of the eigenvalue
#' spectrum and the "signal" part.
#'
#' @section Usage:
#' cutoff(p) %::% matrix : numeric
#'
#' cutoff(p, es, estimator) %::% matrix : list : RandomMatrixFilter : numeric
#'
#' @name cutoff
#' @aliases RandomMatrixFilter MaximumLikelihoodFit
#' @param p A random matrix
#' @return The eigenvalue associated with the upper bound of the noise
#' spectrum
#' @examples
#' cutoff(rmatrix(WignerModel(100)))
cutoff(p) %::% matrix : numeric
cutoff(p) %as% {
  es <- eigen(p, symmetric=TRUE)
  cutoff(p, es, RandomMatrixFilter())
}

# Provide a default implementation
cutoff(p, es, estimator) %::% matrix : eigen : RandomMatrixFilter : numeric
cutoff(p, es, estimator) %when% {
  ! estimator %hasa% fit.fn
} %as% {
  fitter <- MaximumLikelihoodFit(hint=c(1,1))
  estimator$fit.fn <- function(es) fit.density(es, fitter)
  cutoff(p, es, estimator)
}

cutoff(p, es, estimator) %::% matrix : eigen : RandomMatrixFilter : numeric
cutoff(p, es, estimator) %as% {
  fit <- estimator$fit.fn(es)
  Q <- fit$par[1]
  lambda.1 <- es$values[1]
  mp.var <- 1 - lambda.1/length(es$values)
  lambda.plus <- qmp(1, svr=Q, var=mp.var)
  lambda.plus
}

# This is for backwards compatibility since newer versions of R use an
# 'eigen' class for the result of eigen instead of 'list'
cutoff(p, es, estimator) %::% matrix : list : RandomMatrixFilter : numeric
cutoff(p, es, estimator) %as% {
  class(es) <- c('eigen', class(es))
  cutoff(p, es, estimator)
}

#' Fit the eigenvalue spectrum to model
#'
#' This only works for fitting the Marcenko-Pastur distribution. It's
#' also not designed for external use.
#'
#' @name fit.density
#' @param lambda Eigenvalues
#' @param fitter A fit function
#' @return Optimal parameters
#' @examples
#' \dontrun{
#' m <- rmatrix(WishartModel(50, 200))
#' es <- eigen(m)
#' fit.density(es, MaximumLikelihoodFit(hint=c(1,1)))
#' }
fit.density(lambda, fitter) %::% eigen : MaximumLikelihoodFit : list
fit.density(lambda, fitter) %as% {
  really.big <- 100000000000
  x <- lambda$values
  fn <- function(theta)
  {
    s <- sum(-log(dmp(x, svr=theta[1], var=theta[2])))
    s <- ifelse(s == Inf, really.big, s)
    flog.debug("Likelihood value is %s", s)
    s
  }
  optim(fitter$hint, fn)
}

fit.density(lambda, fitter) %::% list : MaximumLikelihoodFit : list
fit.density(lambda, fitter) %as% {
  # For backwards compatibility
  class(lambda) <- 'eigen'
  fit.density(lambda, fitter)
}

