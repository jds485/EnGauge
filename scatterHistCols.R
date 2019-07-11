#Modified scatterHist function from the psych package.
#This one plots kernel smoothed histograms by color.
scatterHist_mod = function (x, y = NULL, smooth = TRUE, ab = FALSE, correl = TRUE, 
          density = TRUE, ellipse = TRUE, digits = 2, method, cex.cor = 1, 
          title = "Scatter plot + histograms", xlab = NULL, ylab = NULL, 
          smoother = FALSE, nrpoints = 0, xlab.hist = NULL, ylab.hist = NULL, 
          grid = FALSE, xlim = NULL, ylim = NULL, x.breaks = 11, y.breaks = 11, 
          x.space = 0, y.space = 0, freq = TRUE, x.axes = TRUE, y.axes = TRUE, 
          size = c(1, 2), col, ...) 
{
  old.par <- par(no.readonly = TRUE)
  n.obs <- sum(!is.na(x))
  if (missing(xlab)) {
    if (!is.null(colnames(x))) {
      xlab = colnames(x)[1]
      ylab = colnames(x)[2]
    }
    else {
      xlab = "V1"
      ylab = "V2"
    }
  }
  if (is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  }else {
    if (!is.null(dim(x))) {
      x <- x[, 1, drop = TRUE]
      if (!is.null(colnames(y))) 
        ylab <- colnames(y)
      if (!is.null(dim(y))) {
        y <- y[, 1, drop = TRUE]
      }
    }
  }
  xrange <- range(x, na.rm = TRUE)
  yrange <- range(y, na.rm = TRUE)
  if (missing(xlim)) 
    xlim <- xrange
  if (missing(ylim)) 
    ylim <- yrange
  x.breaks <- seq(xlim[1], xlim[2], (xlim[2] - xlim[1])/x.breaks)
  y.breaks <- seq(ylim[1], ylim[2], (ylim[2] - ylim[1])/y.breaks)
  xhist <- hist(x, breaks = x.breaks, plot = FALSE)
  yhist <- hist(y, breaks = y.breaks, plot = FALSE)
  nf <- layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), 
               c(3, 1), c(1, 3), TRUE)
  par(mar = c(5, 4, 1, 1))
  if (smoother) {
    smoothScatter(x, y, nrpoints = nrpoints, xlim = xlim, 
                  ylim = ylim, xlab = xlab, ylab = ylab, )
  }else {
    plot(x, y, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col,
         ...)
  }
  if (grid) 
    grid()
  if (ab) 
    abline(lm(y ~ x))
  if (smooth) {
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok]), col = "red")
  }
  if (ellipse) {
    ellipses(x, y, add = TRUE, size = size)
  }
  par(mar = c(0.75, 4, 2, 1))
  #Add other colors to the histogram, if requested
  lst = list()
  ymaxs = 0
  for (co in 1:length(unique(col))){
    al = density(x[col == unique(col)[co]], na.rm = TRUE)
    al$x = ((al$x - min(xhist$breaks))/(max(xhist$breaks) - min(xhist$breaks)))*length(x.breaks)
    lst = c(lst, list(al))
    ymaxs = max(ymaxs, al$y, na.rm=TRUE)
  }
  if (freq) {
    mp <- barplot(xhist$counts, axes = x.axes, space = x.space, 
                  xlab = xlab.hist)
  }else {
    mp <- barplot(xhist$density, axes = x.axes, space = x.space, 
                  xlab = xlab.hist, col = NA, ylim = c(0, ymaxs))
  }
  tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", adjust = 1.2), 
              silent = TRUE)
  if (class(tryd) != "try-error") {
    d$x <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(xhist$breaks))/(max(xhist$breaks) - 
                                                                       min(xhist$breaks))
    if (freq) 
      d$y <- d$y * max(xhist$counts/xhist$density, na.rm = TRUE)
    if (density) 
      lines(d)
  }
  #Add other colors to the histogram, if requested
  for (co in 1:length(unique(col))){
    polygon(lst[[co]], col = adjustcolor(unique(col)[co], alpha.f = 0.3), density = -1, lwd = 2, border = unique(col)[co])
  }
  title(title)
  
  par(mar = c(5, 0.5, 1, 2))
  #Add other colors to the histogram, if requested
  lst = list()
  ymaxs = 0
  for (co in 1:length(unique(col))){
    al = density(y[col == unique(col)[co]], na.rm = TRUE)
    temp = al$y
    al$y = ((al$x - min(yhist$breaks))/(max(yhist$breaks) - min(yhist$breaks)))*length(y.breaks)
    al$x = temp
    lst = c(lst, list(al))
    ymaxs = max(ymaxs, al$x, na.rm=TRUE)
  }
  if (freq) {
    mp <- barplot(yhist$counts, axes = y.axes, space = y.space, 
                  horiz = TRUE, ylab = ylab.hist)
  }else {
    mp <- barplot(yhist$density, axes = y.axes, space = y.space, 
                  horiz = TRUE, ylab = ylab.hist, col = NA, xlim = c(0,ymaxs))
  }
  tryd <- try(d <- density(y, na.rm = TRUE, bw = "nrd", adjust = 1.2), 
              silent = TRUE)
  if (class(tryd) != "try-error") {
    temp <- d$y
    d$y <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(yhist$breaks))/(max(yhist$breaks) - 
                                                                       min(yhist$breaks))
    d$x <- temp
    if (freq) 
      d$x <- d$x * max(yhist$counts/yhist$density, na.rm = TRUE)
    if (density) 
      lines(d)
  }
  for (co in 1:length(unique(col))){
    polygon(lst[[co]], col = adjustcolor(unique(col)[co], alpha.f = 0.3), density = -1, lwd = 2, border = unique(col)[co])
  }
  par(mar = c(1, 1, 1, 1))
  if (correl) {
    plot(1, 1, type = "n", axes = FALSE)
    med.x <- median(x, na.rm = TRUE)
    med.y <- median(y, na.rm = TRUE)
    if (missing(method)) 
      method <- "pearson"
    r = (cor(x, y, use = "pairwise", method = method))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    if (missing(cex.cor)) {
      cex <- 0.75/strwidth(txt)
    }
    else {
      cex <- cex.cor
    }
    text(1, 1, txt, cex = cex)
  }
  par(old.par)
}