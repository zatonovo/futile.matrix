# TODO: remove any ids in columns/rows that don't appear in the triplet form as
# this breaks the CSR construction algorithm.
read.matrix <- function(date, name, format, ..., filter.fn=NULL)
{
  date <- td(date)
  client <- re.options('client')
  if (is.null(client)) { stop("re.option client is required") }
  home <- re.options('home')
  if (is.null(home)) { home = "." }

  path <- sprintf(format, home, client,fd(date))
  zz <- gzfile(path)

  # First row contains column names
  logger.debug('Loading ids')
  all.ids <- strsplit(readLines(zz, n=2), ',', fixed=TRUE)
  col.ids <- all.ids[[1]]
  row.ids <- all.ids[[2]]

  #col.ids <- scan(zz, what='character', sep=',', nlines=1, quiet=TRUE)
  col.ids <- col.ids[col.ids != '']

  # Second row contains row names
  #logger.debug('Loading row ids')
  #row.ids <- scan(zz, what='character', sep=',', nlines=1, skip=1, quiet=TRUE)
  row.ids <- row.ids[row.ids != '']

  logger.debug('Reading triplet representation of %s', name)
  zz <- gzfile(path)
  pts <- read.csv(zz, colClasses=c('character','character','numeric'), skip=2,
    header=FALSE)
  colnames(pts) <- c('row.id','col.id','value')
  logger.debug('Got raw %s: [%s,%s]', name, nrow(pts), ncol(pts))
  if (!is.null(filter.fn))
  {
    logger.info("Applying filter to raw data")
    pts <- filter.fn(pts)
    row.ids <- filter.fn(row.ids)
    logger.debug('Filtered %s: [%s,%s]', name, nrow(pts), ncol(pts))
  }
  #tryCatch(close(zz), finally=logger.debug("Closed file"))

  # This guarantees that the computed ias are monotonically increasing
  row.ids <- row.ids[order(row.ids)]
  col.ids <- col.ids[order(col.ids)]
  pts <- pts[order(pts$row.id,pts$col.id),]

  # Remove records in id lists that are not present in the actual data body.
  # This is necessary to conform to CSR construction rules.
  row.ids <- row.ids[row.ids %in% pts$row.id]
  col.ids <- col.ids[col.ids %in% pts$col.id]

  logger.debug("Output matrix will be [%s,%s]", length(row.ids),length(col.ids))

  m <- matrix.assign.sparse(pts, row.ids, col.ids, ...)
  logger.debug('Assigned values to %s', name)
  logger.debug('Converting to non-sparse matrix')
  m <- as.matrix(m)
  logger.debug('Attaching names')
  colnames(m) <- toupper(col.ids)
  rownames(m) <- toupper(row.ids)
  logger.debug('Done with matrix')
  m
}

# Creates a sparse matrix in CSR format based on a triplet input
# Note that inputs must all be ordered otherwise this will violate the CSR
# spec during construction
matrix.assign.sparse <- function(source, row.ids, col.ids,
  row.block=500000, col.block=500000)
{
  require('SparseM')
  msg.col <- "Calculating column indexes for %s elements (map size: %s)"
  msg.row <- "Calculating row indexes for %s elements (map size: %s)"
  msg.lookup <- "Looking up block %s/%s [%s:%s] of source"

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
    logger.debug(msg.lookup, idx, num.pieces, inf, sup)
    map[source[inf:sup]]
  }

  logger.debug(msg.col, nrow(source), m)
  col.map <- as.integer(1:length(col.ids))
  names(col.map) <- col.ids
  col.idx.list <- apply(array(1:num.pieces), 1, lookup, col.map, source$col.id)
  if (! class(col.idx.list) %in% 'list') col.idx.list <- list(col.idx.list)
  # This is the map of (s)ids to their corresponding column number in the
  # column definition. Values are the index position and the names are the
  # (s)ids. The ordering needs to be consistent with the row.idx
  col.idx <- do.call(c, col.idx.list)

  logger.debug(msg.row, nrow(source), n)
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

  logger.debug("Calculating ias", n)
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
      logger.warn("Skipping bad row match at index %s", idx)
      next
    }
    if (r == prev) next
    if (idx <= 0)
    {
      logger.warn("Bad index value encountered: %s", idx)
    }
    idxes[r] <- idx
    prev <- r
  }, finally=logger.debug("prev=%s, r=%s",prev,r))
  idxes[length(idxes)] <- nrow(source) + 1

  logger.debug("Creating sparse matrix with dimensions [%s,%s]", n,m)
  m.ra <- source$value
  m.ja <- as.integer(col.idx)
  m.ia <- as.integer(idxes)
  d <- as.integer(c(n,m))

  new('matrix.csr', ra=m.ra, ja=m.ja, ia=m.ia, dimension=d)
}

matrix.assign.full <- function(source, row.ids, col.ids)
{
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
  logger.debug('Creating raw rows for matrix')
  all.rows <- lapply(row.ids, fn, source)
  logger.debug('Binding %s rows to create matrix', length(all.rows))
  m <- matrix(do.call(rbind, all.rows), nrow=length(all.rows), byrow=TRUE)
  rownames(m) <- row.ids
  colnames(m) <- col.ids
  m
}

