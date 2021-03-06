#Modified flow-duration curve function from hydroTSM package
fdc.default = function (x, lQ.thr = 0.7, hQ.thr = 0.2, plot = TRUE, log = "y", 
          main = "Flow Duration Curve", xlab = "% Time flow equalled or exceeded", 
          ylab = "Q, [m3/s]", ylim = NULL, yat = c(0.01, 0.1, 1), 
          xat = c(0.01, 0.025, 0.05), col = "black", pch = 1, lwd = 1, 
          lty = 1, cex = 0.4, cex.axis = 1.2, cex.lab = 1.2, leg.txt = NULL, 
          leg.cex = 1, leg.pos = "topright", verbose = TRUE, thr.shw = TRUE, 
          new = TRUE, ...) 
{
  Qposition <- function(x, Q) {
    Q.dist <- abs(x - Q)
    Q.index <- which.min(Q.dist)
    return(Q.index)
  }
  x <- as.numeric(x)
  x.old <- x
  x <- sort(x)
  x.zero.index <- which(x == 0)
  nzeros <- length(x.zero.index)
  ind <- match(x.old, x)
  n <- length(x)
  dc <- rep(NA, n)
  dc[1:n] <- sapply(1:n, function(j, y) {
    dc[j] <- length(which(y >= y[j]))
  }, y = x)
  dc <- dc/n
  if (plot) {
    dc.plot <- dc
    if (log == "y") {
      if (nzeros > 0) {
        x <- x[-x.zero.index]
        dc.plot <- dc.plot[-x.zero.index]
        if (verbose) 
          message("[Note: all 'x' equal to zero (", 
                  nzeros, ") will not be plotted ]")
      }
    }
    if (is.null(ylim)) 
      ylim <- range(x, na.rm = TRUE)
    if (((log == "y") | (log == "xy") | (log == "yx")) & 
        min(ylim) == 0) {
      tmp <- x
      tmp[which(tmp == 0)] <- NA
      ylim[1] <- min(tmp, na.rm = TRUE)
    }
    if (new) {
      plot(dc.plot, x, xaxt = "n", yaxt = "n", type = "o", 
           col = col, pch = pch, lwd = lwd, lty = lty, 
           cex = cex, cex.axis = cex.axis, cex.lab = cex.lab, 
           main = main, xlab = xlab, ylab = ylab, ylim = ylim, 
           log = log, ...)
    }
    else lines(dc.plot, x, xaxt = "n", type = "o", col = col, 
               pch = pch, lwd = lwd, lty = lty, cex = cex)
    ylabels <- pretty(ylim)
    if ((log == "y") | (log == "xy") | (log == "yx")) {
      #MODIFIED HERE
      #Used to be union(ylabels, yat)
      ylabels <- yat
    }
    axis(side = 2, at = ylabels, cex.axis = cex.axis, labels = ylabels)
    xpos <- seq(0, 1, by = 0.05)
    xlabels <- seq(0, 1, by = 0.1)
    if ((log == "x") | (log == "xy") | (log == "yx")) {
      xpos <- union(xat, xpos)
      xlabels <- union(xat, xlabels)
    }
    axis(side = 1, at = xpos, cex.axis = cex.axis, labels = FALSE)
    axis(side = 1, at = xlabels, cex.axis = cex.axis, labels = paste(100 * 
                                                                       xlabels, "%", sep = ""))
    if (!is.na(lQ.thr)) 
      abline(v = lQ.thr, col = "grey", lty = 3, lwd = 2)
    if (!is.na(hQ.thr)) 
      abline(v = hQ.thr, col = "grey", lty = 3, lwd = 2)
    if (!is.null(leg.txt)) 
      legend(x = leg.pos, legend = leg.txt, cex = leg.cex, 
             col = col, pch = pch, lwd = lwd, lty = lty, 
             bty = "n")
    if (thr.shw) {
      x.lQ <- x[Qposition(dc.plot, lQ.thr)]
      x.hQ <- x[Qposition(dc.plot, hQ.thr)]
      legend("bottomleft", c(paste("Qhigh.thr=", round(x.hQ, 
                                                       2), sep = ""), paste("Qlow.thr=", round(x.lQ, 
                                                                                               2), sep = "")), cex = 0.8, bty = "n")
    }
  }
  dc <- dc[ind]
  return(dc)
}