#' Plot simulated polynomial hazards.
#'
#' \code{ggpoly} uses ggplot2 to plot simulated relative hazards from a simpoly class object.
#' @param obj a simpoly object
#'
#'
#' @returns a ggplot2 object.
#'
#' @import ggplot2
#' @export

ggpoly <- function(obj, xlab = NULL, ylab = NULL, title = NULL, xbreaks = NULL, xlabels = NULL, smoother = "auto", colour = "#A6CEE3", spalette = "Set1", leg.name = "", lsize = 2, psize = 1, palpha = 0.1, ...){
  if (!inherits(obj, "simpoly")){
  	stop("must be a simpoly object")
  }

  objdf <- data.frame(obj$X, obj$RH)
  names(objdf) <- c("X", "RH")
  ggplot(objdf, aes(X, RH)) +
    geom_point(shape = 21, alpha = I(palpha), size = psize) +
    geom_smooth(method = smoother, size = lsize, se = FALSE) +
    geom_hline(aes(yintercept = 1), linetype = "dotted") +
    scale_colour_brewer(palette = spalette, name = leg.name) +
    scale_y_continuous()+
    scale_x_continuous() +
    xlab(xlab) + ylab(ylab) +
    ggtitle(title) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw(base_size = 15)

}