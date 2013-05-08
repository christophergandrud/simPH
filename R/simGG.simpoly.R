#' Plot simulated polynomial quantities of interest from Cox Proportional Hazards Models.
#'
#' \code{simGG.simpoly} uses ggplot2 to plot simulated relative quantities of interest from a \code{simpoly} class object.
#' @param obj a simpoly object.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient
#' @param spalette colour palette for use in \code{qi = "Hazard Rate"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "First Difference"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz +
#'        hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum + 
#'        natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + 
#'        condavg3 + orderent + stafcder, data = CarpenterFdaData)
#' 
#' # Simulate simpoly First Difference
#' Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference", 
#'            pow = 3, Xj = seq(1, 150, by = 5))
#'
#' # Simulate simpoly Hazard Ratio with spin probibility interval
#' Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio", 
#'            pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)
#' 
#' # Plot simulations
#' simGG(Sim1)
#' simGG(Sim2)
#'
#' @seealso \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
#'
#' @return a ggplot2 object.
#'
#' @import ggplot2
#' @method simGG simpoly
#' @S3method simGG simpoly

simGG.simpoly <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, ...)
{
  if (!inherits(obj, "simpoly")){
      stop("must be a simpoly object")
    }
    # Find quantity of interest
    qi <- class(obj)[[2]]

    # Create horizontal line (NULL value)
    HL <- 1
    if (qi == "First Difference"){
      HL <- 0
    }    
    # Create y-axis label
    if (is.null(ylab)){
      ylab <- paste(qi, "\n")
    } else {
      ylab <- ylab
    }

    # Subset simlinear object & create a data frame of important variables
  if (qi == "Hazard Rate"){
    colour <- NULL
    if (is.null(obj$strata)){
      objdf <- data.frame(obj$time, obj$QI, obj$HRValue)
      names(objdf) <- c("Time", "HRate", "HRValue")
    } else if (!is.null(obj$strata)) {
    objdf <- data.frame(obj$time, obj$QI, obj$strata, obj$HRValue)
    names(objdf) <- c("Time", "HRate", "Strata", "HRValue")
    }
    if (!is.null(from)){
      objdf <- subset(objdf, Time >= from)
      }
      if (!is.null(to)){
        objdf <- subset(objdf, Time <= to)
      }
  } else if (qi == "Hazard Ratio" | qi == "Relative Hazard" | qi == "First Difference"){
      objdf <- data.frame(obj$Xj, obj$QI)
      names(objdf) <- c("Xj", "QI")
  }

  # Plot
    if (qi == "Hazard Rate"){
      if (!is.null(obj$strata)) {
        ggplot(objdf, aes(x = Time, y = HRate, colour = factor(HRValue))) +
          geom_point(alpha = I(palpha), size = psize) +
          geom_smooth(method = smoother, size = lsize, se = FALSE) +
          facet_grid(.~ Strata) +
          xlab(xlab) + ylab(ylab) +
          scale_colour_brewer(palette = spalette, name = leg.name) +
          ggtitle(title) +
          guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
      } else if (is.null(obj$strata)){
          ggplot(objdf, aes(Time, HRate, colour = factor(HRValue))) +
            geom_point(shape = 21, alpha = I(palpha), size = psize) +
            geom_smooth(method = smoother, size = lsize, se = FALSE) +
            scale_colour_brewer(palette = spalette, name = leg.name) +
            xlab(xlab) + ylab(ylab) +
            ggtitle(title) +
            guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            theme_bw(base_size = 15)
    }
  } else if (qi == "Hazard Ratio" | qi == "Relative Hazard" | qi == "First Difference"){
      ggplot(objdf, aes(Xj, QI)) +
          geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
          geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
          geom_hline(aes(yintercept = HL), linetype = "dotted") +
          xlab(xlab) + ylab(ylab) +
          ggtitle(title) +
          guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
  } 
}