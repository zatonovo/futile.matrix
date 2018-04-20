#'  A collection of matrix manipulation functions
#'
#' This package provides functions for working with random matrices.
#' It also provides various convenience functions for examining data within 
#' matrices as well as some optimized functions for reading matrices in 
#' various formats.
#'
#' \tabular{ll}{
#' Package: \tab futile.matrix\cr
#' Type: \tab Package\cr
#' Version: \tab 1.2.7\cr
#' Date: \tab 2018-04-20\cr
#' License: \tab LGPL-3\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Random matrix ensembles can be created using this package. It's also
#' possible to fit the Marcenko-Pastur distribution to Wishart matrices,
#' enabling you to isolate the noise portion of the eigenvalue spectrum.
#' 
#' @name futile.matrix-package
#' @aliases futile.matrix-package futile.matrix
#' @docType package
#' @exportPattern "^[^\\.]"
#' @import utils lambda.r lambda.tools futile.logger RMTstat
#' @author Brian Lee Yung Rowe <r@zatonovo.com> 
#' @seealso  \code{\link{select}}, \code{\link{expand}}, \code{\link{read.matrix}} 
#' @references 
#' The Distribution Functions of Random Matrix Theory, Craig A. Tracy, UC Davis
#' \url{http://www.math.ucsc.edu/research/rmtg.html}
#'
#' Introduction to the Random Matrix Theory: Gaussian Unitary Ensemble and Beyond
#' \url{http://arxiv.org/abs/math-ph/0412017v2}
#'
#' Tyler's M-Estimator, Random Matrix Theory, and Generalized Elliptical 
#' Distributions with Applications to Finance
#' \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=1287683}
#'
#' @keywords  package attribute logic
#' @examples
#' # Generate a random ensemble
#' m <- rmatrix(WishartModel(100,400))
#' 
#' # Select sub-matrices
#' library(datasets)
#' select(swiss, "Rive")
#' select(swiss, col.pat='^E')
#' select(swiss, "Rive", '^E') <- -1

#' dimnames <- list( c(rownames(swiss), 'Zermat', 'Zurich', 'Geneva'),
#'  c(colnames(swiss), 'Age','Hair.Color') )
#' my.swiss <- expand(swiss, dimnames)
NULL
