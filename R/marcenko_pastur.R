# Get the cutoff for the noise portion of the spectrum, which is essentially
# lambda+. Depending ont the technique, there are different ways of achieving
# this.
cutoff <- function()
{
}

# Calculate and plot the theoretical density distribution
# e.values - The eigenvalues to plot the density against. This can really be any
#   point on the xaxis.
mp.theory <- function(Q, sigma, e.values=NULL, steps=200)
{
  # Plot a range of values
  if (is.null(e.values)) { e.values <- mp.lambdas(Q,sigma,steps) }
  rho <- mp.rho(Q,sigma, e.values)

  if (anylength(e.values) > 1)
  {
    l.min <- mp.eigen.min(Q,sigma)
    l.max <- mp.eigen.max(Q,sigma)
    xs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
    main <- paste('Marcenko-Pastur Distribution for Q',Q,'and sigma',sigma)
    plot(xs,rho, xlim=c(0,6), type='l', main=main)
  }
  rho
}


# Generate eigenvalues for theoretical Marcenko-Pastur distribution
mp.lambdas <- function(Q,sigma, steps)
{
  l.min <- mp.eigen.min(Q,sigma)
  l.max <- mp.eigen.max(Q,sigma)

  logger(DEBUG, sprintf("min eigenvalue:%s",l.min))
  logger(DEBUG, sprintf("max eigenvalue:%s",l.max))

  evs <- seq(round(l.min-1), round(l.max+1), (l.max-l.min)/steps)
  evs[evs < l.min] <- l.min
  evs[evs > l.max] <- l.max

  evs
}

# This provides the theoretical density for a set of eigenvalues. These are
# really just points along the x axis for which the eigenvalue density is
# desired.
# e.values can be a vector of eigen values or a single eigen value
mp.rho <- function(Q,sigma, e.values)
{
  l.min <- mp.eigen.min(Q,sigma)
  l.max <- mp.eigen.max(Q,sigma)

  k <- Q / (2*pi*sigma^2)
  rho <- k * sqrt(pmax(0, (l.max-e.values)*(e.values-l.min)) ) / e.values
  rho
}

# Get maximum eigenvalue specified by Marcenko-Pastur
mp.eigen.max <- function(Q,sigma) { sigma^2 * (1 + sqrt(1/Q))^2 }

# Get minimum eigenvalue specified by Marcenko-Pastur
mp.eigen.min <- function(Q,sigma) { sigma^2 * (1 - sqrt(1/Q))^2 }

