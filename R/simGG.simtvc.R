#' Plot simulated time-varying hazard ratios or stratified time-varying hazard rates from a simtvc class object using ggplot2
#' 
#' \code{simGG.simtvc} uses ggplot2 to plot the simulated hazards from a \code{simtvc} class object using ggplot2. 
#' @param obj a simtvc class object
#' @param from numeric time to start the plot from.
#' @param to numeric time to plot to.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param title the plot's main title
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient
#' @param spalette colour palette for use in \code{qi = "Hazard Rate"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "Relative Hazard"} or \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "Relative Hazard"} or \code{qi = "First Difference"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#'
#' @return a ggplot2 object
#' @details Plots either a time varying hazard ratio or the hazard rates for multiple strata. Currently the strata legend labels need to be changed manually (see \code{\link{revalue}} in the \link{plyr} package) in the \code{simtvc} object with the \code{strata} component. Also, currently the x-axis tick marks and break labels must be adjusted manually for non-linear functions of time. 
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no effect, for time-varying hazard ratio graphs. No line is created for hazard rates.
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#' 
#' # Load survival package
#' library(survival)
#' 
#' # Create natural log time interactions
#' Golubtvc <- function(x){
#'   assign(paste0("l", x), tvc(GolubEUPData, b = x, tvar = "end", tfun = "log"))
#' }
#' 
#' GolubEUPData$Lcoop <-Golubtvc("coop")
#' GolubEUPData$Lqmv <- Golubtvc("qmv")
#' GolubEUPData$Lbacklog <- Golubtvc("backlog")
#' GolubEUPData$Lcodec <- Golubtvc("codec")
#' GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
#' GolubEUPData$Lthatcher <- Golubtvc("thatcher") 
#' 
#' # Run Cox PH Model
#' M1 <- coxph(Surv(begin, end, event) ~ 
#'             qmv + qmvpostsea + qmvpostteu + 
#'             coop + codec + eu9 + eu10 + eu12 +
#'             eu15 + thatcher + agenda + backlog +
#'             Lqmv + Lqmvpostsea + Lcoop + Lcodec +
#'             Lthatcher + Lbacklog, 
#'          data = GolubEUPData,
#'          ties = "efron")
#'          
#' # Create simtvc object for Relative Hazard
#' Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
#'                    tfun = "log", from = 80, to = 2000, 
#'                    Xj = 1, by = 15, ci = 0.99)
#' 
#' # Create simtvc object for First Difference  
#' Sim2 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
#'                   qi = "First Difference", 
#'                   tfun = "log", from = 80, to = 2000, 
#'                   Xj = 1, by = 15, ci = 0.99, spin = TRUE)
#' 
#' # Create simtvc object for Hazard Ratio  
#' Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
#'                   qi = "Hazard Ratio", Xj = c(191, 229), 
#'                   Xl = c(0, 0),
#'                   tfun = "log", from = 80, to = 2000, 
#'                   by = 15, ci = 0.99)
#'                   
#' # Create plots
#' simGG(Sim1)
#' simGG(Sim2)
#' simGG(Sim3, leg.name = "Comparision", from = 1200)

#' @import ggplot2
#' @export
#' @method simGG simtvc
#' @S3method simGG simtvc
#'
#' @references Licht, Amanda A. 2011. ''Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.'' Political Analysis 19: 227-43.

simGG.simtvc <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1)
{
  if (!inherits(obj, "simtvc")){
    stop("must be a simtvc object")
  }

  # Find quantity of interest
  qi <- class(obj)[[2]]

  # Create y-axis label
  if (is.null(ylab)){
    ylab <- paste(qi, "\n")
  } else {
    ylab <- ylab
  }

  # Subset simtvc object & create data frame of important variables
  if (qi == "Hazard Rate"){
    colour <- NULL
    if (is.null(obj$strata)){
      objdf <- data.frame(obj$RealTime, obj$QI, obj$HRValue)
      names(objdf) <- c("Time", "HRate", "HRValue")
    } else if (!is.null(obj$strata)) {
      objdf <- data.frame(obj$RealTime, obj$QI, obj$strata, obj$HRValue)
      names(objdf) <- c("Time", "HRate", "Strata", "HRValue")
    }
  } else if (qi == "Hazard Ratio"){
      objdf <- data.frame(obj$RealTime, obj$QI, obj$Comparison)
      names(objdf) <- c("Time", "QI", "Comparison")
  } else if (qi == "Relative Hazard"){
      objdf <- data.frame(obj$RealTime, obj$QI, obj$Xj)
      names(objdf) <- c("Time", "QI", "Xj")
  } else if (qi == "First Difference"){
      objdf <- data.frame(obj$RealTime, obj$QI, obj$Comparison)
      names(objdf) <- c("Time", "QI", "Comparison")
  }

  # Keep certain times
  if (!is.null(from)){
    objdf <- subset(objdf, Time >= from)
  }
  if (!is.null(to)){
    objdf <- subset(objdf, Time <= to)
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
  } else if (qi == "Hazard Ratio"){
      ggplot(objdf, aes(x = Time, y = QI, colour = factor(Comparison))) +
        geom_point(alpha = I(palpha), size = psize) +
        geom_smooth(method = smoother, size = lsize, se = FALSE) +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        scale_colour_brewer(palette = spalette, name = leg.name) +
        ggtitle(title) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw(base_size = 15)
  } else if (qi == "Relative Hazard"){
      ggplot(objdf, aes(x = Time, y = QI, colour = factor(Xj))) +
        geom_point(alpha = I(palpha), size = psize) +
        geom_smooth(method = smoother, size = lsize, se = FALSE) +
        geom_hline(aes(yintercept = 1), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        scale_colour_brewer(palette = spalette, name = leg.name) +
        ggtitle(title) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw(base_size = 15)
  } else if (qi == "First Difference"){
      ggplot(objdf, aes(Time, QI, group = Comparison)) +
        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
        geom_hline(aes(yintercept = 0), linetype = "dotted") +
        xlab(xlab) + ylab(ylab) +
        scale_colour_brewer(palette = spalette, name = leg.name) +
        ggtitle(title) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme_bw(base_size = 15)
  }
}