#' Plot simulated linear multiplicative interactions from Cox Proportional Hazards Models
#'
#' \code{simGG.siminteract} uses \link{ggplot2} to plot the quantities of interest from \code{siminteract} objects, including marginal effects, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a \code{siminteract} class object
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param from numeric time to start the plot from. Only relevant if \code{qi = "Hazard Rate"}.
#' @param to numeric time to plot to. Only relevant if \code{qi = "Hazard Rate"}.
#' @param title the plot's main title.
#' @param smoother what type of smoothing line to use to summarize the plotted coefficient.
#' @param spalette colour palette. Not relevant for \code{qi = "Marginal Effect"}. Default palette is \code{"Set1"}. See \code{\link{scale_colour_brewer}}.
#' @param leg.name name of the stratified hazard rates legend. Only relevant if \code{qi = "Hazard Rate"}.
#' @param lcolour character string colour of the smoothing line. The default is hexadecimal colour \code{lcolour = '#2B8CBE'}. Only relevant if \code{qi = "Marginal Effect"}.
#' @param lsize size of the smoothing line. Default is 2. See \code{\link{ggplot2}}.
#' @param pcolour character string colour of the simulated points. Default is hexadecimal colour \code{pcolour = '#A6CEE3'}. Only relevant if \code{qi = "Marginal Effect"}.
#' @param psize size of the plotted simulation points. Default is \code{psize = 1}. See \code{\link{ggplot2}}.
#' @param palpha point alpha (e.g. transparency). Default is \code{palpha = 0.05}. See \code{\link{ggplot2}}.
#' @param ribbons logical specifies whether or not to use summary ribbons of the simulations rather than plotting every simulation value as a point. If \code{lines = TRUE} a plot will be created with shaded areas ('ribbons') for the minimum and maximum simulation values (i.e. the middle interval set with \code{qi} in \code{\link{coxsimInteract}}) as well as the central 50% of this area. It also plots a line for the median value of the full area, so values in \code{smoother} are ignored. One of the key advantages of using ribbons rather than points is that it creates plots with smaller file sizes.
#' @param ... Additional arguments. (Currently ignored.)
#'
#' @return a \code{gg} \code{ggplot} class object
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
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx", 
#'							X2 = seq(2, 115, by = 2), nsim = 100)
#'
#' # Plot quantities of interest
#' simGG(Sim1)
#'
#' ## dontrun
#' # Change the order of the covariates to make a more easily
#' # interpretable hazard ratio graph.
#' # M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal, data = CarpenterFdaData)
#'
#' # Simulate Hazard Ratio of lethal for multiple values of prevgenx
#' # Sim2 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal", 
#' #                    X1 = seq(2, 115, by = 2),
#' #                    X2 = c(0, 1),
#' #                    qi = "Hazard Ratio", ci = 0.9)
#'                     
#' # Simulate First Difference
#' # Sim3 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal", 
#' #                       X1 = seq(2, 115, by = 2),
#' #                       X2 = c(0, 1),
#' #                       qi = "First Difference", spin = TRUE)
#'                        
#' # Plot quantities of interest
#' # simGG(Sim1, xlab = "\nprevgenx", ylab = "Marginal Effect of lethal\n")
#' # simGG(Sim2)
#' # simGG(Sim3)
#'
#' @details Uses \link{ggplot2} to plot the quantities of interest from \code{siminteract} objects, including marginal effects, first differences, hazard ratios, and hazard rates. If there are multiple strata, the quantities of interest will be plotted in a grid by strata.
#' Note: A dotted line is created at y = 1 (0 for first difference), i.e. no effect, for time-varying hazard ratio graphs. No line is created for hazard rates.
#'
#'
#' Note: if \code{qi = "Hazard Ratio"} or \code{qi = "First Difference"} then you need to have choosen more than one fitted value for \code{X1} in \code{\link{coxsimInteract}}. 
#'
#' @import ggplot2
#' @method simGG siminteract
#' @S3method simGG siminteract
#' 
#' @seealso \code{\link{coxsimInteract}}, \code{\link{simGG.simlinear}}, and \code{\link{ggplot2}}
#' @references Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. ''Understanding Interaction Models: Improving Empirical Analyses.'' Political Analysis 14(1): 63-82.
#'
#' Keele, Luke. 2010. ''Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.'' Political Analysis 18(2): 189-205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.'' American Journal of Political Science 46(3): 490-505.

simGG.siminteract <- function(obj, from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, ribbons = FALSE, ...)
{
	Time <- QI <- HRValue <- X1 <- X2 <- NULL
	if (!inherits(obj, "siminteract")){
    	stop("must be a siminteract object")
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
		      ggplot(obj, aes(x = Time, y = QI, colour = factor(HRValue))) +
		        geom_point(alpha = I(palpha), size = psize) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE) +
		        facet_grid(.~ Strata) +
		        xlab(xlab) + ylab(ylab) +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
	    	} else if (is.null(obj$strata)){
		      	ggplot(obj, aes(Time, QI, colour = factor(HRValue))) +
		        	geom_point(shape = 21, alpha = I(palpha), size = psize) +
			        geom_smooth(method = smoother, size = lsize, se = FALSE) +
			        scale_colour_brewer(palette = spalette, name = leg.name) +
			        xlab(xlab) + ylab(ylab) +
			        ggtitle(title) +
			        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			        theme_bw(base_size = 15)
			}
		} 
		else if (qi == "Marginal Effect"){
			ggplot(obj, aes(X2, QI)) +
			    geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
		        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +  
		        geom_hline(aes(yintercept = 0), linetype = "dotted") + 
			    xlab(xlab) + ylab(ylab) +
			    ggtitle(title) +
			    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			    theme_bw(base_size = 15)
		} 
		else if (qi == "First Difference"){
			X1Unique <- obj[!duplicated(obj[, "X1"]), ]
			if (nrow(X1Unique) <= 1){
				message("X1 must have more than one fitted value.")
			} else {
				ggplot(obj, aes(X1, QI, colour = factor(X2), group = factor(X2))) +
			        geom_point(shape = 21, alpha = I(palpha), size = psize) +
			        geom_smooth(method = smoother, size = lsize, se = FALSE) +
			        geom_hline(aes(yintercept = 0), linetype = "dotted") +
			        scale_colour_brewer(palette = spalette, name = leg.name) +
			        xlab(xlab) + ylab(ylab) +
			        ggtitle(title) +
			        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			        theme_bw(base_size = 15)
		    }
		} 
		else if (qi == "Hazard Ratio"){
			X1Unique <- obj[!duplicated(obj[, "X1"]), ]
			if (nrow(X1Unique) <= 1){
				message("X1 must have more than one fitted value.")
			} else {
				ggplot(obj, aes(X1, QI, colour = factor(X2), group = factor(X2))) +
			        geom_point(shape = 21, alpha = I(palpha), size = psize) +
			        geom_smooth(method = smoother, size = lsize, se = FALSE) +
			        geom_hline(aes(yintercept = 1), linetype = "dotted") +
			        scale_colour_brewer(palette = spalette, name = leg.name) +
			        xlab(xlab) + ylab(ylab) +
			        ggtitle(title) +
			        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			        theme_bw(base_size = 15)
		    }
	    }
	}
	# Plot points ribbons
	else if (isTRUE(ribbons)){
		suppressWarnings(
		if (qi == "Hazard Rate"){
		  	if (!is.null(obj$strata)) {
			obj <- MinMaxLines(df = obj, hr = TRUE, strata = TRUE)
			ggplot(obj, aes(x = Time, y = HRate, colour = factor(HRValue), fill = factor(HRValue))) +
				geom_line(size = lsize, alpha = palpha) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha, linetype = 0) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, linetype = 0) +
				facet_grid(. ~ Strata) +
				xlab(xlab) + ylab(ylab) +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        scale_fill_brewer(palette = spalette, name = leg.name) +
				ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
			theme_bw(base_size = 15)
    	} else if (is.null(obj$Strata)){
			obj <- MinMaxLines(df = obj, hr = TRUE)
	      	ggplot(obj, aes(Time, Median, colour = factor(HRValue), fill = factor(HRValue))) +
		        geom_line(size = lsize) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), ailpha = palpha, linetype = 0) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, linetype = 0) +
		        scale_colour_brewer(palette = spalette, name = leg.name) +
		        scale_fill_brewer(palette = spalette, name = leg.name) +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
			}
		} 
		else if (qi == "Marginal Effect"){
			obj <- MinMaxLines(df = obj, byVars = c("X2"))
			ggplot(obj, aes(X2, Median)) +
		        geom_line(size = lsize, colour = lcolour) +
				geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha, fill = pcolour) +
				geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, fill = pcolour) +
	        	geom_hline(aes(yintercept = 1), linetype = "dotted") +
		        xlab(xlab) + ylab(ylab) +
		        ggtitle(title) +
		        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		        theme_bw(base_size = 15)
		} 
		else if (qi == "Hazard Ratio" | qi == "First Difference"){
			X1Unique <- obj[!duplicated(obj[, "X1"]), ]
			if (nrow(X1Unique) <= 1){
				message("X1 must have more than one fitted value.")
			} else {
				obj <- MinMaxLines(df = obj, byVars = c("X1", "X2"))
				ggplot(obj, aes(X1, Median, colour = factor(X2), fill = factor(X2))) +
			        geom_line(size = lsize) +
					geom_ribbon(aes(ymin = Lower50, ymax = Upper50), alpha = palpha, linetype = 0) +
					geom_ribbon(aes(ymin = Min, ymax = Max), alpha = palpha, linetype = 0) +
					geom_hline(aes(yintercept = 1), linetype = "dotted") +
					scale_colour_brewer(palette = spalette, name = leg.name) +
			        scale_fill_brewer(palette = spalette, name = leg.name) +
			        xlab(xlab) + ylab(ylab) +
			        ggtitle(title) +
			        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
					theme_bw(base_size = 15)
		    }
	    }
	    )
	}
}