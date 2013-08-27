#' Plot simulated linear time-constant quantities of interest from Cox Proportional Hazards Models
#'
#' \code{simGG.simlinear} uses \link{ggplot2} to plot the quantities of interest from \code{simlinear} objects, including relative hazards, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a \code{simlinear} class object.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title.
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient.
#' @param spalette colour palette for use in \code{qi = "Hazard Rate"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "First Difference"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "First Difference"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#' @param ribbons logical specifies whether or not to use summary ribbons of the simulations rather than plotting every simulation value as a point. If \code{lines = TRUE} a plot will be created with shaded areas ('ribbons') for the minimum and maximum simulation values (i.e. the middle interval set with \code{qi} in \code{\link{coxsimLinear}}) as well as the central 50% of this area. It also plots a line for the median value of the full area, so values in \code{smoother} are ignored. One of the key advantages of using ribbons rather than points is that it creates plots with smaller file sizes.
#' @param ... Additional arguments. (Currently ignored.)
#'
#' @return a \code{gg} \code{ggplot} class object
#'
#' @examples
#'
#' ## dontrun
#' # Load survival package
#' # library(survival)
#' # Load Carpenter (2002) data
#' # data("CarpenterFdaData")
#' 
#' # Estimate survival model
#' # M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
#' #            deathrt1 + acutediz + hosp01  + hhosleng +
#' #            mandiz01 + femdiz01 + peddiz01 + orphdum +
#' #            vandavg3 + wpnoavg3 + condavg3 + orderent +
#' #            stafcder, data = CarpenterFdaData)
#' 
#' # Simulate and plot Hazard Ratios for stafcder variable
#' # Sim1 <- coxsimLinear(M1, b = "stafcder", 
#' #						qi = "Hazard Ratio", 
#' #						Xj = seq(1237, 1600, by = 2), spin = TRUE)
#' # simGG(Sim1)
#' 
#' ## dontrun
#' # Simulate and plot Hazard Rate for stafcder variable
#' # Sim2 <- coxsimLinear(M1, b = "stafcder", nsim = 100,
#' #						qi = "Hazard Rate", Xj = c(1237, 1600))
#' # simGG(Sim2)
#'
#' @details Uses \link{ggplot2} to plot the quantities of interest from \code{simlinear} objects, including relative hazards, first differences, hazard ratios, and hazard rates. If there are multiple strata, the quantities of interest will be plotted in a grid by strata.
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no effect, for time-varying hazard ratio graphs. No line is created for hazard rates.
#'
#' @import ggplot2
#' @method simGG simlinear
#' @S3method simGG simlinear
#'
#' @seealso \code{\link{coxsimLinear}}, \code{\link{simGG.simtvc}}, and \code{\link{ggplot2}}
#' @references Licht, Amanda A. 2011. ''Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.'' Political Analysis 19: 227-43.
#'
#' Keele, Luke. 2010. ''Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.'' Political Analysis 18(2): 189-205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.'' American Journal of Political Science 46(3): 490-505.

simGG.simlinear <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, ribbons = FALSE, ...)
{
	Time <- HRate <- HRValue <- Xj <- QI <- NULL
	if (!inherits(obj, "simlinear")){
    	stop("must be a simlinear object")
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
    # Constrict time period to plot for hazard rate
    if (qi == "Hazard Rate"){   
	    if (!is.null(from)){
			obj <- subset(obj, Time >= from)
		}
		if (!is.null(to)){
	        	obj <- subset(obj, Time <= to)
	    }
    }

	# Plot points
	if (!isTRUE(ribbons)){
		if (qi == "Hazard Rate"){
	  	if (!is.null(obj$strata)) {
			ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue))) +
				geom_point(alpha = I(palpha), size = psize) +
				geom_smooth(method = smoother, size = lsize, se = FALSE) +
				facet_grid(.~ Strata) +
				xlab(xlab) + ylab(ylab) +
				scale_colour_brewer(palette = spalette, name = leg.name) +
				ggtitle(title) +
				guides(colour = guide_legend(override.aes = list(alpha = 1))) +
				theme_bw(base_size = 15)
    	} else if (is.null(obj$strata)){
	      	ggplot(obj, aes(Time, HRate, colour = factor(HRValue))) +
	        	geom_point(shape = 21, alpha = I(palpha), size = psize) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE) +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
			}
		} else if (qi == "First Difference"){
			ggplot(obj, aes(Xj, QI)) +
		        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
		        geom_hline(aes(yintercept = 0), linetype = "dotted") +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
		} else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
		ggplot(obj, aes(Xj, QI)) +
	        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
		}
	}
	# Plot ribbons
	else if (isTRUE(ribbons)){
		############ Incomplete ############
		if (qi == "Hazard Rate"){
	  	if (!is.null(obj$strata)) {
			obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
			ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue))) +
				geom_point(alpha = I(palpha), size = psize) +
				geom_smooth(method = smoother, size = lsize, se = FALSE) +
				facet_grid(.~ Strata) +
				xlab(xlab) + ylab(ylab) +
				scale_colour_brewer(palette = spalette, name = leg.name) +
				ggtitle(title) +
				guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			theme_bw(base_size = 15)
    	} else if (is.null(obj$strata)){
			obj <- MinMaxLines(df = obj, hr = TRUE)
	      	ggplot(obj, aes(Time, Median, colour = factor(HRValue), fill = factor(HRValue))) +
		        geom_line(size = lsize, alpha = I(palpha)) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha) +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
		}
		} else if (qi == "First Difference"){
			obj <- MinMaxLines(df = obj)
			ggplot(obj, aes(Xj, Median)) +
		        geom_line(size = lsize, alpha = I(palpha), colour = lcolour) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha, fill = pcolour) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, fill = pcolour) +
		        geom_hline(aes(yintercept = 0), linetype = "dotted") +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
		} else if (qi == "Hazard Ratio" | qi == "Relative Hazard"){
			obj <- MinMaxLines(df = obj)
			ggplot(obj, aes(Xj, Median)) +
		        geom_line(size = lsize, colour = lcolour) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha, fill = pcolour) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, fill = pcolour) +
	        	geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
		}
	}
}