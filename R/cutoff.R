

#' Filter type for Wigner matrices
#'
#' Uses semi-circle to fit
#' 
#' @references
#' [1] https://terrytao.wordpress.com/tag/wigner-semi-circular-law/
#' @name WignerFilter
#' @param fit.fn The fit function
WignerFilter(fit.fn) %as% list(fit.fn=fit.fn)
WignerFilter() %as% {

}


#' Filter type for Wishart matrices
#'
#' Uses Marcenko-Pastur distribution for the fit.
#' 
#' @name WishartFilter
#' @param fit.fn The fit function
WishartFilter(fit.fn) %as% list(fit.fn=fit.fn)
WishartFilter() %as% {
  fn <- function(x) sum(-log(dmp(x, svr=theta[1], var=theta[2])))
  fitter <- MaximumLikelihoodFit(like.fn=fn, hint=c(1,1))
  fit.fn <- function(es) fit.density(es, fitter)
  WishartFilter(fit.fn)
}


#' Obsolete
#'
#' @name RandomMatrixFilter
#' @param \dots Arguments
RandomMatrixFilter(...) %as% WishartFilter(...)

#' Type that represents MLE fit
#'
#' @name MaximumLikelihoodFit
#' @param like.fn The likelihood function
#' @param hint A hint for the fit function
MaximumLikelihoodFit(like.fn, hint) %as% list(like.fn=like.fn, hint=hint)


#' Find the upper cutoff of the noise part of the spectrum
#'
#' Get the cutoff for the noise portion of the spectrum, which is essentially
#' lambda+ for Wishart matrices. Depending on the technique, there are 
#' different ways of achieving this.
#'
#' @name cutoff
#' @param p A random matrix
#' @examples
#' # This is for demonstration only. Normally you would pass an empirical
#' # correlation matrix
#' cutoff(rmatrix(WishartFilter(100)))
cutoff(p) %::% matrix : numeric
cutoff(p) %as% {
  filter <- .get_filter(p)
  if (is.null(filter))
    stop("Unable to infer matrix type. Please specify a filter explicitly.")
  cutoff(p, filter())
}


cutoff(p, estimator) %::% matrix : list : WishartFilter : numeric
cutoff(p, estimator) %as% {
  es <- eigen(p, symmetric=TRUE)
  fit <- estimator$fit.fn(es)
  Q <- fit$par[1]
  lambda.1 <- es$values[1]
  mp.var <- 1 - lambda.1/length(es$values)
  lambda.plus <- qmp(1, svr=Q, var=mp.var)
  lambda.plus
}


cutoff(p, estimator) %::% matrix : list : WignerFilter : numeric
cutoff(p, estimator) %as% {
  es <- eigen(p)
  fit <- estimator$fit.fn(es)
  Q <- fit$par[1]
  lambda.1 <- es$values[1]
  mp.var <- 1 - lambda.1/length(es$values)
  lambda.plus <- qmp(1, svr=Q, var=mp.var)
  lambda.plus
}


#' Estimate parameters for eigenvalue spectrum
#'
#' Normally this is called within a filter and not called directly.
#'
#' NOTE: When evaluating the fit of estimated parameters, don't use the KS 
#' test [1]. Instead, use [2, 3].
#'
#' 
#' @references
#' [1] http://astrostatistics.psu.edu/su07/R/html/stats/html/ks.test.html
#' [2] http://www-cdf.fnal.gov/physics/statistics/notes/cdf1285_KS_test_after_fit.pdf
#' [3] http://journals.ametsoc.org/doi/pdf/10.1175/MWR3326.1
#' 
#' @examples
#' model <- WishartModel(50, 200)
#' mat <- rmatrix(model)
#' es <- eigen(mat)
#' fit.density(es, MaximumLikelihoodFit(hint=c(1,1)))
fit.density(lambda, fitter) %::% list : MaximumLikelihoodFit : list
fit.density(lambda, fitter) %as% {
  really.big <- 100000000000
  x <- lambda$values
  fn <- function(theta)
  {
    s <- fitter$like.fn(x)
    s <- ifelse(s == Inf, really.big, s)
    flog.debug("Likelihood value is %s", s)
    s
  }
  optim(fitter$hint, fn)
}

