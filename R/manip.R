
#' Peek inside a matrix or vector
#'
#' Peek is a simple utility to conveniently look at a portion of a matrix.
#' This is similar to head and tail but provides a 2-dimensional slice 
#' instead of a complete row. This is useful for debugging and inspecting
#' matrices.
#'
#' @param x Any object that supports subsetting
#' @param upper The upper bound in the subsetting
#' @param lower The lower bound in the subsetting
#' @return A subset of the original matrix, data.frame, etc.
#' @keywords array
#' @examples
#' m <- matrix(c(1,3,4,2, 5,10,11,2, 3,42,8,22, 23,15,3,8), ncol=4)
#' peek(m, 2)
peek <- function(x, upper=5, lower=1)
{
  if (is.null(dim(x)))
  {
    my.upper <- min(upper, anylength(x))
    return(x[lower:my.upper])
  }

  upper.row <- min(upper, anylength(x))
  upper.col <- min(upper, ncol(x))
  return(x[lower:upper.row,lower:upper.col])
}


#' Expand a matrix to larger dimensions, filling in new entries
#'
#' Provides a convenient mechanism for expanding the dimensions of a matrix 
#' with a specified default value. This is particularly useful if the 
#' matrix needs to match dimensions with another matrix. Expand can take 
#' either another matrix or a set of dimnames to grow the matrix.
#' In either case, the expanded matrix will have dimensions that match the 
#' explicitly or implicitly specified dimnames.
#'
#' To properly expand \code{m} to target, the rownames and colnames of 
#' \code{m} must be a strict subset of target's rownames and colnames.
#' If this requirement is not
#' satisfied, the behavior is currently undefined (although most likely an 
#' error will result). In the future, the behavior might be configurable to 
#' drop those rows/columns that are not in target's rownames/colnames.
#'
#' In general, \code{expand} tries to err on the side of accomodation,
#' although the implementation is incomplete. If target is a list, then 
#' the format is the same as when constructing a matrix and passing 
#' dimnames as an argument. Currently, only a list or a matrix are supported.
#' If a list, target[[0]] represent the row names and target[[1]] are 
#' the column names.
#' This could be relaxed in the future to any object that has a rownames 
#' and colnames. 
#'
#' Note that a current limitation/feature in expand is that it orders the
#' resulting matrix by rows and columns. More precise control needs to be
#' provided here, with the default being the ordering of the rows and columns
#' conforming to target.
#'
#' TODO: Consistency check to ensure all rownames/colnames of m are a subset
#' of target
#'
#' @param m A matrix to expand
#' @param target A list containing the dimnames or a matrix that contains 
#'  dimnames
#' @param default The default value to use for the new entries
#' @return The expanded matrix with rows and columns corresponding to the 
#'  target dimnames
#' @keywords array
#' @examples
#' rows.m <- c('row.1', 'row.2', 'row.3')
#' cols.m <- c('col.1', 'col.2')
#' rows.n <- c(rows.m, 'row.4')
#' cols.n <- c(cols.m, 'col.3')
#' m <- matrix(c(1,4,7,2,5,8), ncol=2, dimnames=list(rows.m,cols.m))
#' n <- matrix(c(1,4,7,10,2,5,8,11,3,6,9,12), ncol=3,
#'   dimnames=list(rows.n,cols.n))
#' expand(m, n)
expand <- function(m, target, default=0) {
  if ('list' %in% class(target)) { dims <- target }
  else if ('matrix' %in% class(target))
  {
    dims <- list(rownames(target), colnames(target))
  }
  else { stop("Argument target must be either a matrix or a list") }

  # Build out rows to new dimensions
  filler.row <- matrix(
    rep(default, ncol(m) * (anylength(dims[[1]]) - nrow(m))), ncol=ncol(m),
    dimnames=list(setdiff(dims[[1]], rownames(m)), colnames(m)) )

  full <- rbind(m, filler.row)
  # Add all other columns
  filler.col <- matrix(
    rep(default, nrow(full) * (anylength(dims[[2]]) - ncol(full))),
    nrow=nrow(full),
    dimnames=list(rownames(full), setdiff(dims[[2]], colnames(full))) )
  full <- cbind(full, filler.col)
  # Sort
  arrange(full)
}

#' Select a portion of a matrix based on a regular expression of the row and/or
#' column names.
#'
#' Extract a subset of a matrix based on regex patterns on either
#' the rownames, the colnames or both. Once this subset has been selected,
#' assignments can be made following standard consistency rules.
#'
#' Oftentimes it is useful to get at a specific subset of data within a matrix.
#' In large matrices, it can be cumbersome to access specific rows and/or
#' columns using traditional subsetting methods, particularly if it is a complex
#' set that is to be extracted. \code{select} provides regex searching on named
#' matrices to access portions of a matrix that satisfy the regex. Note that
#' \code{select} will work for data.frames as well.
#'
#' It is possible to assign values to the selected subset as a means to modify
#' the original matrix. Standard consistency rules must be satisfied for any
#' assignment operations.
#'
#' @param m A matrix from which to select a subset 
#' @param row.pat A regular expression to use for rownames
#' @param col.pat A regular expression to use for colnames
#' @param ... Additional arguments to pass to grep
#' @return  A matrix containing all rows and columns that satisfy the 
#'  patterns given. If no values match, then an empty matrix will be returned.
#' @aliases select<-
#' @examples
#' library(datasets)
#' select(swiss, "Rive")
#'
#' select(swiss, col.pat="E", fixed=TRUE)
#'
#' select(swiss, row.pat='^[A-T]', col.pat="^E")
#'
#' select(swiss, "Rive", "Ed") <- min(select(swiss, "^[^R]", "Ed"))
select <- function(m, row.pat=NULL, col.pat=NULL, ...) {
  out <- m
  if (! is.null(row.pat)) {
    out <- out[grep(row.pat, rownames(out), ...), , drop=FALSE]
  }

  if (! is.null(col.pat)) {
    out <- out[, grep(col.pat, colnames(out), ...), drop=FALSE]
  }
  out
}


# Select a portion of a matrix based on a regular expression and assign the
# subset to value. Dimensional integrity is required, otherwise an error will
# result.
"select<-" <- function(m, row.pat=NULL, col.pat=NULL, ..., value) {
  if (is.null(row.pat) & is.null(col.pat))
  { stop("Either row.pat or col.pat must be set") }

  if (! is.null(row.pat) & is.null(col.pat)) {
    m[grep(row.pat, rownames(m), ...),] <- value
  }
  else if (is.null(row.pat) & ! is.null(col.pat)) {
    m[, grep(col.pat, colnames(m), ...)] <- value
  }
  else {
    rows <- grep(row.pat, rownames(m), ...)
    cols <- grep(col.pat, colnames(m), ...)
    m[rows,cols] <- value
  }
  invisible(m)
}


#' Order matrix elements based on rownames and/or colnames
#'
#' This is a convenience function for ordering rows and columns within a
#' matrix.
#'
#' To ensure proper operations are performed, the ordering of rows and columns
#' within matrix operands must be consistent. Arrange conveniently performs
#' this ordering.
#'
#' In the future, a comparator will be added so that custom orderings can be
#' applied to the function.
#'
#' @param m A matrix whose rows and/or columns are to be arranged 
#' @param order.rows Whether rows are to be ordered. Defaults to TRUE 
#' @param order.cols Whether columns are to be ordered. Defaults to TRUE 
#' @param comparator A function to define the ordering. Currently unused 
#' @return A matrix whose rows and/or columns have been ordered. By default,
#'  both rows and columns are ordered.
#' @examples
#' library(datasets)
#' arrange(state.x77)
arrange <- function(m, order.rows=TRUE, order.cols=TRUE, comparator=NULL) {
  #if (is.null(col.ids)) col.ids <- colnames(m)
  #if (is.null(row.ids)) row.ids <- rownames(m)

  if (order.rows & order.cols)
    m[order(rownames(m)),order(colnames(m)), drop=FALSE]
  else if(order.cols)
    m[,order(colnames(m)), drop=FALSE]
  else if(order.rows)
    m[order(rownames(m)), , drop=FALSE]
}

