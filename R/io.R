.onLoad <- function(libname, pkgname)
{
  .resetLogger()
}

.resetLogger <- function()
{
  require(futile.logger)
  .log <<- getLogger("futile.matrix")
}

# TODO: remove any ids in columns/rows that don't appear in the triplet form as
# this breaks the CSR construction algorithm.
# Alternate test
# rows <- paste('row',1:10000, sep='.')
# cols <- paste('col',1:10000, sep='.')
# mat <- cbind(sample(rows, 5000), sample(cols, 5000), rnorm(5000))
# write.csv(mat, file='triplet.csv', row.names=FALSE, col.names=FALSE, quote=FALSE)

.fileConnection <- function(file)
{
  if ('connection' %in% class(file)) { zz <- file }
  else if (length(grep("\\.gz(ip)?$", file)) > 0) { zz <- gzfile(file) }
  else { zz <- file(file) }
  zz
}

.processHeader <- function(file, row.idx, col.idx)
{
  zz <- .fileConnection(file)
  .log(DEBUG,'Loading ids')

  count <- 0
  if (! is.null(row.idx)) count <- count + 1
  if (! is.null(col.idx)) count <- count + 1
  all.ids <- strsplit(readLines(zz, n=count), ',', fixed=TRUE)

  idx <- 0
  if (! is.null(row.idx))
  {
    row.ids <- all.ids[[idx]]
    row.ids <- row.ids[row.ids != '']
    idx <- idx + 1
  }
  if (! is.null(col.idx))
  {
    col.ids <- all.ids[[idx]]
    col.ids <- col.ids[col.ids != '']
    idx <- idx + 1
  }

  list(row.ids=row.ids, col.ids=col.ids, skip=idx)
}


# About 20s to load a 10000 x 10000 sparse matrix, 413MB memory
# rows <- paste('row', 1:10000, sep='.')
# cols <- paste('col', 1:10000, sep='.')
# m <- read.matrix('triplet.csv', row.ids=rows, col.ids=cols)
# 
# Compare to slam ~18s, 2.3GB memory
# read.slam <- function(file)
# {
#   m <- read.csv(file, header=FALSE, skip=1)
#   m[,1] <- sub('row.|col.','', m[,1])
#   m[,2] <- sub('row.|col.','', m[,2])
#   m1 <- simple_triplet_matrix(m[,1], m[,2], m[,3], nrow=10000, ncol=10000)
#   as.matrix(m1)
# }
# system.time(m <- read.slam('triplet.csv'))
read.matrix <- function(file, header=FALSE, skip=1, 
  row.ids=NULL, col.ids=NULL,
  colClasses=c('character','character','numeric'), 
  assign.fn=assignMatrixTriplet, filter.fn=NULL, ...)
{
  #logger <- getLogger("futile.matrix")

  if (header)
  {
    as <- .processHeader(file, row.ids, col.ids)
    if (! is.null(as)) attach(as)
  }

  .log(DEBUG,sprintf('Reading triplet representation of %s', file))
  zz <- .fileConnection(file)
  pts <- read.csv(zz, colClasses=colClasses, skip=skip, header=FALSE)
  colnames(pts) <- c('row.id','col.id','value')
  .log(DEBUG,sprintf('Got raw %s: [%s,%s]', file, nrow(pts), ncol(pts)))

  if (is.null(row.ids)) row.ids <- pts[,1]
  if (is.null(col.ids)) col.ids <- pts[,2]

  if (!is.null(filter.fn))
  {
    .log(INFO,"Applying filter to raw data")
    pts <- filter.fn(pts)
    row.ids <- filter.fn(row.ids)
    .log(DEBUG,sprintf('Filtered %s: [%s,%s]', file, nrow(pts), ncol(pts)))
  }
  #tryCatch(close(zz), finally=.log(DEBUG,"Closed file"))

  # This guarantees that the computed ias are monotonically increasing
  row.ids <- row.ids[order(row.ids)]
  col.ids <- col.ids[order(col.ids)]
  pts <- pts[order(pts$row.id,pts$col.id),]

  #m <- assignMatrixSparse(pts, row.ids, col.ids, ...)
  m <- assign.fn(pts, row.ids, col.ids, ...)
  .log(DEBUG,sprintf('Assigned values to %s', file))
  .log(DEBUG,'Converting to dense matrix')
  m <- as.matrix(m)
  .log(DEBUG,'Attaching names')
  colnames(m) <- toupper(col.ids)
  rownames(m) <- toupper(row.ids)
  .log(DEBUG,'Done with matrix')
  m
}

# Creates a sparse matrix in CSR format based on a triplet input
# Note that inputs must all be ordered otherwise this will violate the CSR
# spec during construction
assignMatrixSparse <- function(source, row.ids, col.ids,
  row.block=500000, col.block=500000)
{
  require(futile.logger)
  require(SparseM)
  msg.col <- "Calculating column indexes for %s elements (map size: %s)"
  msg.row <- "Calculating row indexes for %s elements (map size: %s)"
  msg.lookup <- "Looking up block %s/%s [%s:%s] of source"

  # Remove records in id lists that are not present in the actual data body.
  # This is necessary to conform to CSR construction rules.
  row.ids <- row.ids[row.ids %in% source$row.id]
  col.ids <- col.ids[col.ids %in% source$col.id]
  .log(DEBUG,sprintf("Output matrix will be [%s,%s]", length(row.ids),length(col.ids)))

  n <- length(row.ids)
  m <- length(col.ids)

  # This blocking is done for performance reasons. There seems to be some
  # cutoff in vector size in certain operations that causes serious slowdowns.
  block.size <- col.block
  num.pieces <- nrow(source) %/% block.size + 1
  lookup <- function(idx, map, source)
  {
    inf <- block.size * (idx - 1) + 1
    sup <- ifelse(idx == num.pieces, length(source), idx * block.size)
    .log(DEBUG,sprintf(msg.lookup, idx, num.pieces, inf, sup))
    map[source[inf:sup]]
  }

  .log(DEBUG,sprintf(msg.col, nrow(source), m))
  col.map <- as.integer(1:length(col.ids))
  names(col.map) <- col.ids
  col.idx.list <- apply(array(1:num.pieces), 1, lookup, col.map, source$col.id)
  if (! class(col.idx.list) %in% 'list') col.idx.list <- list(col.idx.list)
  # This is the map of (s)ids to their corresponding column number in the
  # column definition. Values are the index position and the names are the
  # (s)ids. The ordering needs to be consistent with the row.idx
  col.idx <- do.call(c, col.idx.list)

  .log(DEBUG,sprintf(msg.row, nrow(source), n))
  block.size <- row.block
  num.pieces <- nrow(source) %/% block.size + 1
  row.map <- as.integer(1:length(row.ids))
  names(row.map) <- row.ids
  row.idx.list <- apply(array(1:num.pieces), 1, lookup, row.map, source$row.id)
  if (! class(row.idx.list) %in% 'list') row.idx.list <- list(row.idx.list)
  # This is the map of (s)ids to their corresponding row number in the row
  # definition. Values are the index position and the names are the (s)ids.
  # The ordering needs to be consistent with the row.idx
  row.idx <- do.call(c, row.idx.list)

  .log(DEBUG,sprintf("Calculating ias", n))
  idxes <- rep(0, n+1)
  idxes[1] <- 1
  idx <- 0
  prev <- 1
  # The idxes become the ias in the CSR format. For this to work, the row ids
  # need to be ordered otherwise the position values are not monotonically
  # increasing which violates the CSR spec.
  tryCatch(
  for (r in row.idx)
  {
    idx <- idx + 1
    if (is.na(r))
    {
      .log(WARN,sprintf("Skipping bad row match at index %s", idx))
      next
    }
    if (r == prev) next
    if (idx <= 0)
    {
      .log(WARN,sprintf("Bad index value encountered: %s", idx))
    }
    idxes[r] <- idx
    prev <- r
  }, finally=.log(DEBUG,sprintf("prev=%s, r=%s",prev,r)))
  idxes[length(idxes)] <- nrow(source) + 1

  .log(DEBUG,sprintf("Creating sparse matrix with dimensions [%s,%s]", n,m))
  m.ra <- source$value
  m.ja <- as.integer(col.idx)
  m.ia <- as.integer(idxes)
  d <- as.integer(c(n,m))

  new('matrix.csr', ra=m.ra, ja=m.ja, ia=m.ia, dimension=d)
}

# Assign the matrix as though it were dense
assignMatrixDense <- function(source, row.ids, col.ids)
{
  require(futile.logger)
  num.cols <- length(col.ids)
  source.rows <- unique(source$row.id)

  fn <- function(row.id, source)
  {
    row <- rep(0, num.cols)
    if (! row.id %in% source.rows) return(row)

    items <- source[source$row.id == row.id,]
    row[col.ids %in% items$col.id] <- items$value
    row
  }
  .log(DEBUG,'Creating raw rows for matrix')
  all.rows <- lapply(row.ids, fn, source)
  .log(DEBUG,sprintf('Binding %s rows to create matrix', length(all.rows)))
  m <- matrix(do.call(rbind, all.rows), nrow=length(all.rows), byrow=TRUE)
  rownames(m) <- row.ids
  colnames(m) <- col.ids
  m
}

# Assign the matrix using the SLAM triplet construction
assignMatrixTriplet <- function(source, row.ids, col.ids, ...)
{
  require(futile.logger)
  require(slam)

  row.idx <- 1:length(row.ids)
  names(row.idx) <- row.ids
  col.idx <- 1:length(col.ids)
  names(col.idx) <- col.ids
 
  i <- row.idx[source$row.id]
  j <- col.idx[source$col.id]
  v <- source$value
  simple_triplet_matrix(i,j,v, nrow=length(row.ids), ncol=length(col.ids), ...)
}


