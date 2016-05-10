#' Perform the conjugate transpose of a matrix
#'
#' Convenience function
#'
#' @section Usage:
#' ct(m) \%::\% matrix : matrix
#'
#' ct(m)
#' 
#' @section Details:
#' This is a convenience function to compute the conjugate transpose. For
#' real-valued matrices, ct(m) = t(m).
#'
#' @name ct
#' @param m A matrix
#' @return THe conjugate transpose of the original matrix
#' @examples
#' x <- matrix(rcomp(16), nrow=4)
#' ct(x)
ct(m) %::% matrix : matrix
ct(m) %as% { Conj(t(m)) }


#' Generate random complex numbers
#'
#' Generate random complex numbers using the specified distribution.
#' By default \code{rnorm} is used.
#'
#' @section Usage:
#' rcomp(n, dist) %::% numeric : Function : complex
#' 
#' rcomp(n, dist=rnorm)
#'
#' @section Details:
#' This function is used primarily to generate random matrices.
#'
#' @name rcomp
#' @param n Length of the output vector
#' @param dist The distribution for the random number genertor
#' @return A vector of random numbers
#' @examples
#' rcomp(10)
#' rcomp(10, runif)
rcomp(n, dist) %::% numeric : Function : complex
rcomp(n, dist=rnorm) %as% {
  complex(real=dist(n), imaginary=dist(n))
}

#' Generation of random matrices
#'
#' Generate various types of random matrices
#'
#' @section Usage:
#' rmatrix(model) %::% WignerModel : matrix
#'
#' rmatrix(model) %::% WishartModel : matrix
#'
#' rmatrix(model) %::% JacobiModel : matrix
#'
#' rmatrix(model)
#'
#' @section Details:
#' Used to generate a random matrix from various families. The idea
#' is to specify a model, which is then used to generate random realizations
#' and also to compute other properties of the matrix.
#'
#' @name rmatrix
#' @aliases symmetric hermitian max_eigen eigenvalues
#' @param model The matrix model to use, which includes the size of
#' the matrix. The model argument must be of type RandomMatrixModel. 
#' Numerous sub-types (e.g. WignerModel, WishartModel) are
#' supported generating the appropriate type of random matrix.
#' @seealso \code{\link{dmatrix}}
#' @examples
#' model <- WignerModel(10)
#' m <- rmatrix(model)
#'
#' \dontrun{
#' e <- Ensemble(20, model)
#' hist(max_eigen(e), freq=FALSE)
#' }

rmatrix(model) %::% WignerModel : matrix
rmatrix(model) %when% {
  !model$real
  # result == ct(result)
} %as% {
  n <- model$n
  x <- matrix(rcomp(n^2), nrow=n) / 2^0.5
  WignerMatrix((x + ct(x)) / sqrt(2 * n), model)
}

# Also known as a Gaussian Unitary Ensemble
rmatrix(model) %::% WignerModel : matrix
rmatrix(model) %when% {
  model$real
  # result == t(result)
} %as% {
  n <- model$n
  x <- matrix(rnorm(n^2), nrow=n)
  WignerMatrix((x + ct(x)) / sqrt(2 * n), model)
}

# For definition of self-dual quaternion and matrix representation,
# http://www.aimath.org/conferences/ntrmt/talks/Mezzadri2.pdf
# http://tonic.physics.sunysb.edu/~verbaarschot/lecture/lecture2.ps

rmatrix(model) %::% WishartModel : matrix
rmatrix(model) %when% {
  !model$real
} %as% {
  n <- model$n
  m <- model$m
  dist.fn <- function(x) rnorm(x, sd=model$sd)
  x <- matrix(rcomp(n * m, dist=dist.fn), nrow=n) / 2^0.5
  WishartMatrix((x %*% ct(x)) / m, model)
}

rmatrix(model) %::% WishartModel : matrix
rmatrix(model) %when% {
  model$real
} %as% {
  n <- model$n
  m <- model$m
  x <- matrix(rnorm(n * m, sd=model$sd), nrow=n)
  WishartMatrix((x %*% t(x)) / m, model)
}


rmatrix(model) %::% JacobiModel : matrix
rmatrix(model) %when% {
  model$real
} %as% {
  n <- model$n
  m1 <- model$m1
  m2 <- model$m2

  x1 <- rmatrix(WishartModel(n,m1, real=model$real))
  x2 <- rmatrix(WishartModel(n,m2, real=model$real))
  JacobiMatrix(solve(x1 + x2) * x1, model)
}


#' Type constructors for random matrices and ensembles of random matrices
#'
#' Provides type constructors for creating random matrices. Various
#' studies can be initiated afterward.
#'
#' @section Usage:
#' RandomMatrixModel(real=TRUE, ...)
#'
#' WignerMatrix(x, model)
#'
#' WishartModel(n, m, sd=1, ...)
#'
#' JacobiModel(n, m1, m2, ...)
#'
#' Ensemble(count, model)
#'
#' @name RandomMatrixModel
#' @aliases WignerModel WignerMatrix WishartModel WishartMatrix JacobiModel JacobiMatrix Ensemble
#' @param real Whether the matrix has real components or not
#' @param n Number of rows
#' @param m Number of columns
#' @param m1 Number of columns
#' @param m2 Number of columns
#' @param sd Standard deviation of the sample population
#' @param count Number of matrices in the ensemble
#' @param model The random matrix model to use
#' @return Returns a model type. Use with \code{\link{rmatrix}} or 
#' \code{\link{Ensemble}} to generate actual matrices.
#' @examples
#' model <- WignerModel(10)
#' m <- rmatrix(model)
#' e <- Ensemble(20, model)
RandomMatrixModel(real=TRUE, ...) %as% list(real=real, ...)

# Random square matrix. Eienvalues form semicircle
WignerModel(n, ...) %as% {
  RandomMatrixModel(n=n, ...)
}

WignerMatrix(x, model) %as% {
  x@n <- model$n
  x
}

# n - variables
# m - observations
# model <- WishartModel(100,500, sd=1)
# hist(eigenvalues(rmatrix(model)))
WishartModel(n, m, sd=1, ...) %as%
{
  RandomMatrixModel(n=n, m=m, Q=m/n, sd=sd, ...)
}

WishartMatrix(x, model) %as% {
  x@n <- model$n
  x@m <- model$m
  x@Q <- model$Q
  x@sd <- model$sd
  x
}

JacobiModel(n, m1, m2, ...) %as%
{
  RandomMatrixModel(n=n, m1=m1, m2=m2, ...)
}

JacobiMatrix(x, model) %as% {
  x@n <- model$n
  x@m1 <- model$m1
  x@m2 <- model$m2
  x
}


Ensemble(count, model) %as%
{
  out <- lapply(seq(count), function(junk) rmatrix(model))
  out@model <- class(model)[1]
  out
}

#' Print a random matrix ensemble
#'
#' Pretty print the ensemble instead of dumping a bunch of matrices
#'
#' @param x An ensemble of random matrices
#' @param ... Reserved for later
#' @return Used for side-effects
#' @seealso \code{\link{Ensemble}}
# @export
#' @S3method print Ensemble
print.Ensemble <- function(x, ...)
{
  cat("\nClass:", attr(x,'model'))
  cat("\nCount:", length(x))
  cat("\nDimensions:", dim(x[[1]]))
  cat("\n")
  invisible(x)
}



eigenvalues(m) %::% matrix : numeric
eigenvalues(m) %as% {
  o <- eigen(m, only.values=TRUE)
  o$values
}

# Example
# en <- Ensemble(50, WignerModel(200))
# hist(max_eigen(en), freq=FALSE)
max_eigen(ensemble) %::% Ensemble : numeric
max_eigen(ensemble) %as% {
  es <- lapply(ensemble, eigen, only.values=TRUE)
  sapply(es, function(x) x$values[1])
}



# Convenience functions
hermitian(n) %as% rmatrix(WignerModel(n=n, real=FALSE))

symmetric(n) %as% rmatrix(WignerModel(n=n, real=TRUE))

