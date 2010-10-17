# References
#   http://arxiv.org/abs/math-ph/0412017v2
#   http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1287683
#   http://web.mit.edu/sea06/agenda/talks/Kuijlaars.pdf
#   http://www.williams.edu/go/math/sjmiller/public_html/BrownClasses/54/handouts/IntroRMT_Math54.pdf
#   http://web.mit.edu/sea06/agenda/talks/Harding.pdf
#   http://dspace.mit.edu/bitstream/handle/1721.1/39670/180190294.pdf
#   http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.cmp/1103842703
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.10.2630
#   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.5.9005
# Resources
#   http://www.math.ucsc.edu/research/rmtg.html
#   The Distribution Functions of Random Matrix Theory, Craig A. Tracy, UC Davis
# Search Terms
#   extreme eigenvalues random matrix feature extraction
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

guard(rmatrix.hermitian, isa(HermitianModel, model))
ensure(rmatrix.hermitian, result == ct(result))
rmatrix.hermitian <- function(rank, model)
{
  m <- diag(rnorm(rank))

  u <- rcomp((rank^2 - rank) / 2)
  m[upper.tri(m)] <- u
  l <- ct(m)
  m[lower.tri(m)] <- l[lower.tri(l)]

  m
}

guard(rmatrix.symmetric, isa(RealSymmetricModel, model))
ensure(rmatrix.symmetric, result == t(result))
rmatrix.symmetric <- function(rank, model)
{
  m <- diag(rnorm(rank))

  u <- rnorm((rank^2 - rank) / 2)
  m[upper.tri(m)] <- u
  l <- t(m)
  m[lower.tri(m)] <- l[lower.tri(l)]

  m
}

# Convenience function
guard(hermitian.default, TRUE)
hermitian.default <- function(rank) rmatrix(rank, create(HermitianModel))

guard(symmetric.default, TRUE)
symmetric.default <- function(rank) rmatrix(rank, create(RealSymmetricModel))


create.Ensemble <- function(T, rank, count, model)
{
  lapply(rep(rank,count), rmatrix, model)
}


# Example
# en <- create(Ensemble, 10, 100, create(HermitianModel))
# hist(max_eigen(en), freq=FALSE)
guard(max_eigen.en, isa(Ensemble,ensemble))
max_eigen.en <- function(ensemble)
{
  es <- lapply(ensemble, eigen, only.values=TRUE)
  sapply(es, function(x) x$values[1])
}
