# References
# The Eigenvalues of Random Matrices Experiments with the Classical Ensembles
#   http://web.mit.edu/18.338/www/handouts/handout3.pdf
# Introduction to the Random Matrix Theory: Gaussian Unitary Ensemble and Beyond
#   http://arxiv.org/abs/math-ph/0412017v2
# Tyler's M-Estimator, Random Matrix Theory, and Generalized Elliptical 
# Distributions with Applications to Finance
#   http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1287683
#   http://web.mit.edu/sea06/agenda/talks/Kuijlaars.pdf
# Distributions of the extreme eigenvalues of the complex Jacobi random matrix
# ensemble
#   http://www.math.washington.edu/~dumitriu/kd_submitted.pdf
#   http://www.williams.edu/go/math/sjmiller/public_html/BrownClasses/54/handouts/IntroRMT_Math54.pdf
#   http://web.mit.edu/sea06/agenda/talks/Harding.pdf
#   http://dspace.mit.edu/bitstream/handle/1721.1/39670/180190294.pdf
#   http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.cmp/1103842703
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.10.2630
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.5.9005
#
# Resources
#   http://www.math.ucsc.edu/research/rmtg.html
#   The Distribution Functions of Random Matrix Theory, Craig A. Tracy, UC Davis
#
# Search Terms
#   extreme eigenvalues random matrix feature extraction
#
# Trading Ideas
#   http://ideas.repec.org/p/lei/ingber/03ai.html
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.7.1041
guard(ct.matrix, is.matrix(m))
ct.matrix <- function(m) Conj(t(m))

guard(rcomp.fun, is.function(dist))
rcomp.fun <- function(n, dist)
{
  sapply(1:n, function(x) complex(real=dist(1), imaginary=dist(1)))
}

guard(rcomp.default, TRUE)
rcomp.default <- function(n) rcomp(n, rnorm)

# Also known as a Gaussian Orthogonal Ensemble
guard(rmatrix.hermitian, isa(WignerModel,model) & !model$real)
ensure(rmatrix.hermitian, result == ct(result))
rmatrix.hermitian <- function(model)
{
  n <- model$n
  x <- matrix(rcomp(n^2), nrow=n) / 2^0.5
  (x + ct(x)) / sqrt(2 * n)
}

# Also known as a Gaussian Unitary Ensemble
guard(rmatrix.symmetric, isa(WignerModel,model) & model$real)
ensure(rmatrix.symmetric, result == t(result))
rmatrix.symmetric <- function(model)
{
  n <- model$n
  x <- matrix(rnorm(n^2), nrow=n)
  (x + ct(x)) / sqrt(2 * n)
}

guard(rmatrix.cwishart, isa(WishartModel,model) & !model$real)
rmatrix.cwishart <- function(model)
{
  n <- model$n
  m <- model$m
  x <- matrix(rcomp(n * m), nrow=n) / 2^0.5
  (x %*% ct(x)) / m
}

guard(rmatrix.rwishart, isa(WishartModel,model) & model$real)
rmatrix.rwishart <- function(model)
{
  n <- model$n
  m <- model$m
  x <- matrix(rnorm(n * m), nrow=n)
  (x %*% t(x)) / m
}


guard(rmatrix.jacobi, isa(JacobiModel,model) & model$real)
rmatrix.jacobi <- function(model)
{
  n <- model$n
  m1 <- model$m1
  m2 <- model$m2

  x1 <- rmatrix(create(WishartModel, n,m1, real=model$real))
  x2 <- rmatrix(create(WishartModel, n,m2, real=model$real))
  solve(x1 + x2) * x1
}

create.RandomMatrixModel <- function(T, real=TRUE) list(real=real)

create.WignerModel <- function(T, n, ...)
{
  o <- create(RandomMatrixModel, ...)
  o$n <- n
  o
}

create.WishartModel <- function(T, n, m, ...)
{
  o <- create(RandomMatrixModel, ...)
  o$n <- n
  o$m <- m
  o
}

create.JacobiModel <- function(T, n, m1, m2, ...)
{
  o <- create(RandomMatrixModel, ...)
  o$n <- n
  o$m1 <- m1
  o$m2 <- m2
  o
}

create.Ensemble <- function(T, rank, count, model)
{
  lapply(rep(rank,count), rmatrix, model)
}


guard(density_lim.wigner_def, isa(WignerModel,model))
density_lim.wigner_def <- function(model) density_lim(-2,2,100, model)

guard(density_lim.wigner, isa(WignerModel,model))
density_lim.wigner <- function(min, max, steps, model)
{
  x <- seq(min, max, length.out=steps)
  sqrt(4 - x^2) / (2 * pi)
}

guard(density_lim.wishart_def, isa(WishartModel,model))
density_lim.wishart_def <- function(model) density_lim(-3,3,100, model)

guard(density_lim.wishart, isa(WishartModel,model))
density_lim.wishart <- function(min, max, steps, model)
{
  x <- seq(min, max, length.out=steps)
  c <- model$n / model$m
  b.neg <- (1 - sqrt(c))^2
  b.pos <- (1 + sqrt(c))^2
  ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  sqrt(ind * (x - b.neg) * (b.pos - x)) / (2*pi*x*c)
}

guard(density_lim.jacobi_def, isa(JacobiModel,model))
density_lim.jacobi_def <- function(model) density_lim(-1,1,100, model)

guard(density_lim.jacobi, isa(JacobiModel,model))
density_lim.jacobi <- function(min, max, steps, model)
{
  lg <- getLogger("futile.matrix")
  lg(WARN, "This function is incomplete")
  x <- seq(min, max, length.out=steps)
  c1 <- model$n / model$m1
  c2 <- model$n / model$m2
  c0 <- c1*x + x^3*c1 - 2*c1*x^2 - c2*x^3 + c2*x^2

  b0 <- c1*x - c2*x - c1 + 2
  b1 <- -2*c2*x^2 + 2*x - 3*c1*x + c1 + c2*x - 1 + 2*c1*x^2
  b2 <- c0

  #num <- (1 - c1)^2
  #b.neg <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 + 2 * sqrt(c1 + c2 - c1 * c2))
  #b.pos <- num / (c1^2 - c1 + 2 + c2 - c1 * c2 - 2 * sqrt(c1 + c2 - c1 * c2))

  #ind <- ifelse(b.neg < x & x < b.pos, 1, 0)
  #sqrt(ind * (x - b.neg) * (x + b.neg)) / (2 * pi * c0)
  sqrt(4*b2*b0 - b1^2) / (2*pi*b2)
}

guard(eigenvalues.matrix, is.matrix(m))
eigenvalues.matrix <- function(m)
{
  o <- eigen(m, only.values=TRUE)
  o$values
}

# Example
# en <- create(Ensemble, 10, 100, create(WignerModel))
# hist(max_eigen(en), freq=FALSE)
guard(max_eigen.en, isa(Ensemble,ensemble))
max_eigen.en <- function(ensemble)
{
  es <- lapply(ensemble, eigen, only.values=TRUE)
  sapply(es, function(x) x$values[1])
}



# Convenience functions
guard(hermitian.default, TRUE)
hermitian.default <- function(rank) rmatrix(rank, create(HermitianModel))

guard(symmetric.default, TRUE)
symmetric.default <- function(rank) rmatrix(rank, create(RealSymmetricModel))

