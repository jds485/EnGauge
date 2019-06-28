#Editing the matplot function to allow for plotting dates on the x-axis instead of the numeric version of the date.
matplotDates = function (x, y, type = "p", lty = 1:5, lwd = 1, lend = par("lend"), 
          pch = NULL, col = 1:6, cex = NULL, bg = NA, xlab = NULL, 
          ylab = NULL, xlim = NULL, ylim = NULL, log = "", ..., add = FALSE, 
          verbose = getOption("verbose")) 
{
  library(zoo)
  paste.ch <- function(chv) paste0("\"", chv, "\"", collapse = " ")
  str2vec <- function(string) {
    if (nchar(string, type = "c")[1L] > 1L) 
      strsplit(string[1L], NULL)[[1L]]
    else string
  }
  xlabel <- if (!missing(x)) 
    deparse(substitute(x))
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  if (missing(x)) {
    if (missing(y)) 
      stop("must specify at least one of 'x' and 'y'")
    else x <- seq_len(NROW(y))
  }
  else if (missing(y)) {
    y <- x
    ylabel <- xlabel
    x <- seq_len(NROW(y))
    xlabel <- ""
  }
  kx <- ncol(x <- as.matrix(x))
  ky <- ncol(y <- as.matrix(y))
  n <- nrow(x)
  if (n != nrow(y)) 
    stop("'x' and 'y' must have same number of rows")
  if (kx > 1L && ky > 1L && kx != ky) 
    stop("'x' and 'y' must have only 1 or the same number of columns")
  if (kx == 1L) 
    x <- matrix(x, nrow = n, ncol = ky)
  if (ky == 1L) 
    y <- matrix(y, nrow = n, ncol = kx)
  k <- max(kx, ky)
  type <- str2vec(type)
  if (is.null(pch)) {
    pch <- c(1L:9L, 0L, letters, LETTERS)
    if (k > length(pch) && any(type %in% c("p", "o", "b"))) 
      warning("default 'pch' is smaller than number of columns and hence recycled")
  }
  else if (is.character(pch)) 
    pch <- str2vec(pch)
  if (verbose) 
    message("matplot: doing ", k, " plots with ", paste0(" col= (", 
                                                         paste.ch(col), ")"), paste0(" pch= (", paste.ch(pch), 
                                                                                     ")"), " ...\n", domain = NA)
  xy <- xy.coords(x, y, xlabel, ylabel, log = log)
  xlab <- if (is.null(xlab)) 
    xy$xlab
  else xlab
  ylab <- if (is.null(ylab)) 
    xy$ylab
  else ylab
  xlim <- if (is.null(xlim)) 
    range(xy$x[is.finite(xy$x)])
  else xlim
  ylim <- if (is.null(ylim)) 
    range(xy$y[is.finite(xy$y)])
  else ylim
  if (length(type) < k) 
    type <- rep_len(type, k)
  if (length(lty) < k) 
    lty <- rep_len(lty, k)
  if (length(lend) < k) 
    lend <- rep_len(lend, k)
  if (length(lwd) < k && !is.null(lwd)) 
    lwd <- rep_len(lwd, k)
  if (length(pch) < k) 
    pch <- rep_len(pch, k)
  if (length(col) < k) 
    col <- rep_len(col, k)
  if (length(bg) < k) 
    bg <- rep_len(bg, k)
  if (is.null(cex)) 
    cex <- 1
  if (length(cex) < k) 
    cex <- rep_len(cex, k)
  ii <- seq_len(k)
  dev.hold()
  on.exit(dev.flush())
  if (!add) {
    ii <- ii[-1L]
    plot(zoo::as.Date(x[, 1L]), y[, 1L], type = type[1L], xlab = xlab, 
         ylab = ylab, xlim = zoo::as.Date(xlim), ylim = ylim, lty = lty[1L], 
         lwd = lwd[1L], lend = lend[1L], pch = pch[1L], col = col[1L], 
         cex = cex[1L], bg = bg[1L], log = log, ...)
  }
  for (i in ii) lines(zoo::as.Date(x[, i]), y[, i], type = type[i], lty = lty[i], 
                      lwd = lwd[i], lend = lend[i], pch = pch[i], col = col[i], 
                      cex = cex[i], bg = bg[i])
  invisible()
}
