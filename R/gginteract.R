#' Plot simulated linear multiplicative interactions.
#'
#' \code{gginteract} uses ggplot2 to plot the quantities of interest from \code{siminteract} objects, including marginal effects, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a simlinear object
#' @param qi character string indicating what quantity of interest you would like to calculate. Can be \code{'Marginal Effect'}, \code{'First Difference'}, \code{'Hazard Ratio'}, or \code{'Hazard Rate'}. Default is \code{qi = 'Marginal Effect'}. 
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient
#' @param spalette colour palette. Not relevant for \code{qi = "Marginal Effect"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "Marginal Effect"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points for relative hazards. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "Marginal Effect"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#' @param ... other arguments passed to specific methods
#' @return a ggplot2 object
#'
#' @examples 
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx, data = CarpenterFdaData)
#' 
#' # Simulate Marginal Effect of lethal for multiple values of prevgenx
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx", X2 = seq(2, 115, by = 2))
#' 
#' # Plot Marginal Effects
#' gginteract(Sim1, qi = "Marginal Effect", xlab = "\nprevgenx", ylab = "Marginal Effect of lethal\n")
#'
#' @description Uses ggplot2 to plot the quantities of interest from \code{siminteract} objects, including marginal effects, first differences, hazard ratios, and hazard rates. If there are multiple strata, the quantities of interest will be plotted in a grid by strata.
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no effect, for time-varying hazard ratio graphs. No line is created for hazard rates.
#'
#'
#' Note: if \code{qi = "Hazard Ratio"} or \code{qi = "First Difference"} then you need to have choosen more than one fitted value for \code{X1} in \code{\link{coxsimInteract}}. 
#'
#' @import ggplot2
#' @export
#' @seealso \code{\link{coxsimInteract}}, \code{\link{gglinear}}, and \code{\link{ggplot2}}
#' @references Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. “Understanding Interaction Models: Improving Empirical Analyses.” Political Analysis 14(1): 63–82.
#'
#' Keele, Luke. 2010. “Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.” Political Analysis 18(2): 189–205.
#'
#' Carpenter, Daniel P. 2002. “Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.” American Journal of Political Science 46(3): 490–505.

gginteract <- function(obj, qi = "Marginal Effect", from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, ...)
{
	if (!inherits(obj, "siminteract")){
    	stop("must be a siminteract object")
    }
    # Create y-axis label
    if (is.null(ylab)){
    	ylab <- paste(qi, "\n")
    } else {
    	ylab <- ylab
    }

    # Subset siminteract object & create a data frame of important variables
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
		colour <- NULL
	  	objdf <- data.frame(obj$X1, obj$X2, obj$HR, obj$Comparison)
	  	names(objdf) <- c("X1", "X2", "HR", "Comparison")
	} else if (qi == "Marginal Effect"){
	  	spalette <- NULL
	  	objdf <- data.frame(obj$X2, obj$HR)
	  	names(objdf) <- c("X2", "HR")
	} else if (qi == "First Difference"){
		colour <- NULL
		objdf <- data.frame(obj$X1, obj$X2, obj$HR)
		names(objdf) <- c("X1", "X2", "FirstDiff")
	}

	# Plot
	if (qi == "Hazard Rate"){
	  	if (!is.null(obj$strata)) {
	      ggplot(objdf, aes(x = Time, y = HRate, colour = factor(HRValue))) +
	        geom_point(alpha = I(palpha), size = psize) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE) +
	        facet_grid(.~ Strata) +
	        scale_y_continuous()+
	        scale_x_continuous() +
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
		        scale_y_continuous()+
		        scale_x_continuous() +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
		}
	} 
	else if (qi == "Marginal Effect"){
		ggplot(objdf, aes(X2, HR)) +
		    geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
		    geom_hline(aes(yintercept = 0), linetype = "dotted") +
		    scale_y_continuous() +
			scale_x_continuous() +    
		    xlab(xlab) + ylab(ylab) +
		    ggtitle(title) +
		    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		    theme_bw(base_size = 15)
	} 
	else if (qi == "First Difference"){
		X1Unique <- objdf[!duplicated(objdf[, "X1"]), ]
		if (nrow(X1Unique) <= 1){
			message("X1 must have more than one fitted value.")
		} else {
			ggplot(objdf, aes(X1, FirstDiff, colour = factor(X2), group = factor(X2))) +
		        geom_point(shape = 21, alpha = I(palpha), size = psize) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE) +
		        geom_hline(aes(yintercept = 0), linetype = "dotted") +
		        scale_y_continuous()+
		        scale_x_continuous() +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
	    }
	} 
	else if (qi == "Hazard Ratio"){
		X1Unique <- objdf[!duplicated(objdf[, "X1"]), ]
		if (nrow(X1Unique) <= 1){
			message("X1 must have more than one fitted value.")
		} else {
			ggplot(objdf, aes(X1, HR, colour = factor(X2), group = factor(X2))) +
		        geom_point(shape = 21, alpha = I(palpha), size = psize) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE) +
		        geom_hline(aes(yintercept = 1), linetype = "dotted") +
		        scale_y_continuous()+
		        scale_x_continuous() +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
	    }
    }
}