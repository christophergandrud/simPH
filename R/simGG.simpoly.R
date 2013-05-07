#' Plot simulated polynomial hazards.
#'
#' \code{simGG.simpoly} uses ggplot2 to plot simulated relative hazards from a \code{simpoly} class object.
#' @param obj a simpoly object
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default is \code{ylab = "Relative Hazard \n"}.
#' @param title the plot's main title
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Works if \code{strata = FALSE}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. 
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#' @param ... other arguments passed to specific methods
#'
#' @examples
#' # Load Carpenter (2002) data
#' # data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' # M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + 
#'   #           acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + 
#'   #           peddiz01 + orphdum + natreg + I(natreg^2) + vandavg3 + 
#'   #           wpnoavg3 + condavg3 + orderent + stafcder, 
#'   #          data = CarpenterFdaData)
#' 
#' # Simulate simpoly class object
#' # Sim1 <- coxsimPoly(M1, b = "natreg", pow = 3, X = seq(1, 150, by = 5))
#' 
#' # Plot simulations
#' # simGG(Sim1)
#'
#' @seealso \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
#'
#' @return a ggplot2 object.
#'
#' @import ggplot2
#' @method simGG simpoly
#' @S3method simGG simpoly

simGG.simpoly <- function(obj, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", lcolour = "#2B8CBE", pcolour = "#A6CEE3",lsize = 2, psize = 1, palpha = 0.1, ...)
{
  if (!inherits(obj, "simpoly")){
  	stop("must be a simpoly object")
  }
    # Find quantity of interest
    qi <- class(obj)[[2]]

    # Create y-axis label
    if (is.null(ylab)){
      ylab <- paste(qi, "\n")
    } else {
      ylab <- ylab
    }

  objdf <- data.frame(obj$Xjl, obj$QI)
  names(objdf) <- c("Xjl", "QI")
  ggplot(objdf, aes(Xjl, QI)) +
        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
        geom_smooth(method = smoother, colour = lcolour, size = lsize, se = FALSE) +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        ggtitle(title) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw(base_size = 15)
}