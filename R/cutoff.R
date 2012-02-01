# Get the cutoff for the noise portion of the spectrum, which is essentially
# lambda+. Depending ont the technique, there are different ways of achieving
# this.

# Assume p is a random matrix
cutoff %when% (p %isa% matrix)
cutoff %as% function(p)
{
  es <- eigen(p, symmetric=TRUE)
  cutoff(p, es, create(RandomMatrixFilter))
}

# Provide a default implementation
cutoff %when% (estimator %isa% RandomMatrixFilter)
cutoff %also% (! estimator %hasa% fit.fn)
cutoff %as% function(p, es, estimator)
{
  fitter <- create(MaximumLikelihoodFit, hint=c(1,1))
  estimator$fit.fn <- function(es) fit.density(es, fitter)
  cutoff(p, es, estimator)
}

cutoff %when% (estimator %isa% RandomMatrixFilter)
cutoff %as% function(p, es, estimator)
{
  fit <- estimator$fit.fn(es)
  Q <- fit$par[1]
  lambda.1 <- es$values[1]
  mp.var <- 1 - lambda.1/length(es$values)
  lambda.plus <- qmp(1, svr=Q, var=mp.var)
  lambda.plus
}


# model <- create(WishartModel, 50, 200)
# mat <- rmatrix(model)
# es <- eigen(mat)
# fit.density(es, create(MaximumLikelihoodFit, hint=c(1,1)))
#
# Don't use ks.test with estimated parameters
#   http://astrostatistics.psu.edu/su07/R/html/stats/html/ks.test.html
# Alternatives:
#   http://www-cdf.fnal.gov/physics/statistics/notes/cdf1285_KS_test_after_fit.pdf
#   http://journals.ametsoc.org/doi/pdf/10.1175/MWR3326.1
# ks.test(es$values, 'pmp', svr=theta[1], var=theta[2])
fit.density %when% (fitter %isa% MaximumLikelihoodFit)
fit.density %as% function(lambda, fitter)
{
  logger <- getLogger('futile.matrix')
  really.big <- 100000000000
  x <- lambda$values
  fn <- function(theta)
  {
    s <- sum(-log(dmp(x, svr=theta[1], var=theta[2])))
    s <- ifelse(s == Inf, really.big, s)
    logger(DEBUG, sprintf("Likelihood value is %s", s))
    s
  }
  optim(fitter$hint, fn)
}

