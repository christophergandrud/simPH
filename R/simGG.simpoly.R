#' Plot simulated polynomial quantities of interest from Cox Proportional Hazards Models
#'
#' \code{simGG.simpoly} uses \link{ggplot2} to plot simulated relative quantities of interest from a \code{simpoly} class object.
#' @param obj a \code{simpoly} class object.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title.
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient.
#' @param spalette colour palette for when there are multiple sets of comparisons to plot. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param legend specifies what type of legend to include (if applicable). The default is \code{legend = "legend"}. To hide the legend use \code{legend = FALSE}. See the \code{\link{discrete_scale}} for more details.
#' @param leg.name name of the legend (if applicable).
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 1. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points or ribbons (when there are not multiple sets of simulations). Default is hexadecimal colour \code{pcolour = '#A6CEE3'}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param alpha point alpha (e.g. transparency) for the points or ribbons. Default is \code{alpha = 0.1}. See \code{\link{ggplot2}}.
#' @param type character string. Specifies how to plot the simulations. Can be \code{points}, \code{lines}, or \code{ribbons}. If points then each simulation value will be plotted. If \code{lines} is chosen then each simulation is plotted using a different line. Note: any simulation with a value along its length that is outside of the specified central interval will be dropped. This is to create a smooth plot. If \code{type = "ribbons"} a plot will be created with shaded areas ('ribbons') for the minimum and maximum simulation values (i.e. the middle interval set with \code{qi} in \code{\link{coxsimSpline}}) as well as the central 50 percent of this area. It also plots a line for the median value of the full area, so values in \code{smoother} are ignored. One of the key advantages of using ribbons rather than points is that it creates plots with smaller file sizes.
#' @param ... Additional arguments. (Currently ignored.)
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + 
#'        deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 + 
#'        femdiz01 + peddiz01 + orphdum + natreg +
#'        I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + 
#'        condavg3 + orderent + stafcder, data = CarpenterFdaData)
#' 
#' # Simulate simpoly First Difference
#' Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference", 
#'            pow = 3, Xj = seq(1, 150, by = 5), nsim = 100)
#'
#' # dontrun
#' # Simulate simpoly Hazard Ratio with spin probibility interval
#' # Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio", 
#' #           pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)
#'
#' # Sim3 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Rate", 
#' #           pow = 3, Xj = c(1, 150))
#' 
#' # Plot simulations
#' simGG(Sim1)
#' # dontrun
#' # simGG(Sim2, ribbons = TRUE)
#' # simGG(Sim3, ribbons = TRUE)
#'
#' @details Uses \link{ggplot2} to plot the quantities of interest from \code{simpoly} objects. 
#'
#' @seealso \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
#'
#' @return a \code{gg} \code{ggplot} class object
#'
#' @import ggplot2
#' @method simGG simpoly
#' @S3method simGG simpoly

simGG.simpoly <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", legend = "legend", leg.name = "", lcolour = "#2B8CBE", lsize = 1, pcolour = "#A6CEE3", psize = 1, alpha = 0.1, type = "points", ...)
{
  Time <- HRValue <- HRate <- Xj <- QI <- Lower50 <- Upper50 <- Min <- Max <- Median <- SimID <- NULL
  if (!inherits(obj, "simpoly")){
    stop("must be a simpoly object")
  }
  if (type == 'ribbons' & smoother != "auto"){
    message("The smoother argument is ignored if ribbons = TRUE. Central tendency summarised with the median.")
  }
  # Find quantity of interest
  qi <- class(obj)[[2]]

  # Create y-axis label
  if (is.null(ylab)){
    ylab <- paste(qi, "\n")
  } else {
    ylab <- ylab
  }
  # Convert obj to data frame
  class(obj) <- "data.frame"

  # Drop simulations that include outliers
  if (type == 'lines'){
    obj <- OutlierDrop(obj)
  }    

  # Constrict time period to plot for hazard rate  
  if (!is.null(from)){
    obj <- subset(obj, Time >= from)
  }
  if (!is.null(to)){
    obj <- subset(obj, Time <= to)
  }

  # Plot points
  if (type == 'points'){
    if (qi == "Hazard Rate"){
      if (!is.null(obj$strata)) {
        ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue))) +
          geom_point(alpha = I(alpha), size = psize) +
          geom_smooth(method = smoother, size = lsize, se = FALSE) +
          facet_grid(.~ Strata) +
          xlab(xlab) + ylab(ylab) +
          scale_colour_brewer(palette = spalette, name = leg.name, guide = legend) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
      } else if (is.null(obj$strata)){
          ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
            geom_point(shape = 21, alpha = I(alpha), size = psize) +
            geom_smooth(method = smoother, size = lsize, se = FALSE) +
            scale_colour_brewer(palette = spalette, name = leg.name, guide = legend) +
            xlab(xlab) + ylab(ylab) +
            ggtitle(title) +
            #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            theme_bw(base_size = 15)
    }
    } else if (qi == "First Difference"){
      ggplot(obj, aes(Xj, QI)) +
          geom_point(shape = 21, alpha = I(alpha), size = psize, colour = pcolour) +
          geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
          geom_hline(aes(yintercept = 0), linetype = "dotted") +
          xlab(xlab) + ylab(ylab) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
    } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
      ggplot(obj, aes(Xj, QI)) +
          geom_point(shape = 21, alpha = I(alpha), size = psize, colour = pcolour) +
          geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
          geom_hline(aes(yintercept = 1), linetype = "dotted") +
          xlab(xlab) + ylab(ylab) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
    }
  }
  # Plot lines
  else if (type == 'lines'){
    if (qi == "Hazard Rate"){
      if (!is.null(obj$strata)) {
        ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue))) +
          geom_line(aes(group = interaction(SimID, factor(HRValue))), alpha = I(alpha), size = psize) +
          geom_smooth(aes(colour = factor(HRValue)), method = smoother, size = lsize, se = FALSE) +
          facet_grid(.~ Strata) +
          xlab(xlab) + ylab(ylab) +
          scale_colour_brewer(palette = spalette, name = leg.name, guide = legend) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
      } else if (is.null(obj$strata)){
          ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
            geom_line(aes(group = interaction(SimID, factor(HRValue))), alpha = I(alpha), size = psize) +
            geom_smooth(aes(colour = factor(HRValue)), method = smoother, size = lsize, se = FALSE) +
            scale_colour_brewer(palette = spalette, name = leg.name, guide = legend) +
            xlab(xlab) + ylab(ylab) +
            ggtitle(title) +
            #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            theme_bw(base_size = 15)
    }
    } else if (qi == "First Difference"){
      ggplot(obj, aes(Xj, QI)) +
          geom_line(aes(group = SimID), alpha = I(alpha), size = psize, colour = pcolour) +
          geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
          geom_hline(aes(yintercept = 0), linetype = "dotted") +
          xlab(xlab) + ylab(ylab) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
    } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
      ggplot(obj, aes(Xj, QI)) +
          geom_line(aes(group = SimID), alpha = I(alpha), size = psize, colour = pcolour) +
          geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
          geom_hline(aes(yintercept = 1), linetype = "dotted") +
          xlab(xlab) + ylab(ylab) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
    }
  }
  # Plot ribbons
  else if (type == 'ribbons'){
    suppressWarnings(
    if (qi == "Hazard Rate"){
      if (!is.null(obj$Strata)) {
        obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
        ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue), fill = factor(HRValue))) +
          geom_line(size = lsize) +
          geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = alpha, linetype = 0) +
          geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha, linetype = 0) +
          facet_grid(. ~ Strata) +
          xlab(xlab) + ylab(ylab) +
          scale_colour_brewer(palette = spalette, name = leg.name, guide = legend) +
          scale_fill_brewer(palette = spalette, name = leg.name, guide = legend) +
          ggtitle(title) +
          #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
          theme_bw(base_size = 15)
      } else if (is.null(obj$Strata)){
      obj <- MinMaxLines(df = obj, hr = TRUE)
          ggplot(obj, aes(Time, Median, colour = factor(HRValue), fill = factor(HRValue))) +
            geom_line(size = lsize) +
            geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = alpha, linetype = 0) +
            geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha, linetype = 0) +
            scale_colour_brewer(palette = spalette, name = leg.name) +
            scale_fill_brewer(palette = spalette, name = leg.name) +
            xlab(xlab) + ylab(ylab) +
            ggtitle(title) +
            #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            theme_bw(base_size = 15)
    }
    } else if (qi == "First Difference"){
      obj <- MinMaxLines(df = obj)
      ggplot(obj, aes(Xj, Median)) +
        geom_line(size = lsize, colour = lcolour) +
        geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = alpha, fill = pcolour) +
        geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha, fill = pcolour) +
        geom_hline(aes(yintercept = 0), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        ggtitle(title) +
        #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            theme_bw(base_size = 15)
    } else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
      obj <- MinMaxLines(df = obj)
      ggplot(obj, aes(Xj, Median)) +
        geom_line(size = lsize, colour = lcolour) +
        geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = alpha, fill = pcolour) +
        geom_ribbon(aes(ymin = Min, ymax = Max), alpha = alpha, fill = pcolour) +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        ggtitle(title) +
        #guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw(base_size = 15)
    }
    )  
  }
}