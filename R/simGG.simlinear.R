#' Plot simulated linear time-constant hazards.
#'
#' \code{simGG.simlinear} uses ggplot2 to plot the quantities of interest from \code{simlinear} objects, including relative hazards, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a simlinear object
#' @param qi character string indicating what quantity of interest you would like to calculate. Can be \code{'Relative Hazard'}, \code{'First Difference'}, \code{'Hazard Ratio'}, or \code{'Hazard Rate'}. Default is \code{qi = 'Relative Hazard'}. 
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient
#' @param spalette colour palette for use in \code{qi = "Hazard Rate"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "Relative Hazard"} or \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "Relative Hazard"} or \code{qi = "First Difference"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#' @param ... other arguments passed to specific methods
#' @return a ggplot2 object
#'
#' @examples
#' # Load survival package
#' library(survival)
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#' 
#' # Estimate survival model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum + natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + condavg3 + orderent + stafcder, data = CarpenterFdaData)
#' 
#' # Simulate and plot Hazard Ratios for stafcder variable
#' Sim1 <- coxsimLinear(M1, b = "stafcder", qi = "Hazard Ratio", Xj = seq(1237, 1600, by = 2))
#' simGG(Sim1, qi = "Hazard Ratio")
#' 
#' # Simulate and plot Hazard Rate for stafcder variable
#' Sim2 <- coxsimLinear(M1, b = "stafcder", qi = "Hazard Rate", Xj = c(1237, 1600))
#' simGG(Sim2, qi = "Hazard Rate")
#'
#' @description Uses ggplot2 to plot the quantities of interest from \code{simlinear} objects, including relative hazards, first differences, hazard ratios, and hazard rates. If there are multiple strata, the quantities of interest will be plotted in a grid by strata.
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no effect, for time-varying hazard ratio graphs. No line is created for hazard rates.
#'
#' @import ggplot2
#' @method simGG simlinear
#' @S3method simGG simlinear
#'
#' @seealso \code{\link{coxsimLinear}}, \code{\link{simGG.simtvc}}, and \code{\link{ggplot2}}
#' @references Licht, Amanda A. 2011. ''Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.'' Political Analysis 19: 227–43.
#' Keele, Luke. 2010. ''Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.'' Political Analysis 18(2): 189–205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.'' American Journal of Political Science 46(3): 490–505.

simGG.simlinear <- function(obj, qi = "Relative Hazard", from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, ...)
{
	if (!inherits(obj, "simlinear")){
    	stop("must be a simlinear object")
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
			objdf <- data.frame(obj$time, obj$HRate, obj$HRValue)
			names(objdf) <- c("Time", "HRate", "HRValue")
		} else if (!is.null(obj$strata)) {
		objdf <- data.frame(obj$time, obj$HRate, obj$strata, obj$HRValue)
		names(objdf) <- c("Time", "HRate", "Strata", "HRValue")
		}
		if (!is.null(from)){
			objdf <- subset(objdf, Time >= from)
  		}
  		if (!is.null(to)){
  			objdf <- subset(objdf, Time <= to)
  		}
	} else if (qi == "Hazard Ratio"){
	  	objdf <- data.frame(obj$Xj, obj$HR)
	  	names(objdf) <- c("Xj", "HR")
	} else if (qi == "Relative Hazard"){
	  	spalette <- NULL
	  	objdf <- data.frame(obj$Xj, obj$HR)
	  	names(objdf) <- c("Xj", "HR")
	} else if (qi == "First Difference"){
		spalette <- NULL
		objdf <- data.frame(obj$Xj, obj$HR, obj$Comparison)
		names(objdf) <- c("Xj", "FirstDiff", "Comparison")
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
	} else if (qi == "Relative Hazard"){
		ggplot(objdf, aes(Xj, HR)) +
		    geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
		    geom_hline(aes(yintercept = 1), linetype = "dotted") +
		    xlab(xlab) + ylab(ylab) +
		    ggtitle(title) +
		    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		    theme_bw(base_size = 15)
	} else if (qi == "First Difference"){
    	ggplot(objdf, aes(Xj, FirstDiff, group = Comparison)) +
        	geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 0), linetype = "dotted") +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
	} else if (qi == "Hazard Ratio"){
		ggplot(objdf, aes(Xj, HR)) +
	        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
    }
}