#' Plot simulated linear time-constant hazards.
#'
#' \code{gglinear} uses ggplot2 to plot the quantities of interest from \code{simlinear} objects.
#'
#' @param obj a simlinear object
#' @param qi character string indicationg what quantity of interest you would like to calculate. Can be \code{'Relative Hazard'}, \code{'First Difference'}, \code{'Hazard Ratio'}, or \code{'Hazard Rate'}. Default is \code{qi = 'Relative Hazard'}. 
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#'
#'
#'
#'
#'
#' @import ggplot2
#' @export
#' @seealso \code{\link{coxsimLinear}} and \code{\link{ggplot2}}
#' @references Licht, Amanda A. 2011. “Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.” Political Analysis 19: 227–43.

gglinear <- function(obj, qi = "Relative Hazard", from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, smoother = "auto", colour = "#A6CEE3", spalette = "Set1", leg.name = "", lsize = 2, psize = 1, palpha = 0.1, ...)
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

    # Subset simlinear object & create data frame of important variables
	if (qi == "Hazard Rate"){
		colour <- NULL
		if (is.null(obj$strata)){
			objdf <- data.frame(obj$time, obj$HRate, obj$Comparison)
			names(objdf) <- c("Time", "HRate", "Comparison")
		} else {
		objdf <- data.frame(obj$time, obj$HRate, obj$strata, obj$Comparison)
		names(objdf) <- c("Time", "HRate", "Strata", "Comparison")
		}
		if (!is.null(from)){
			objdf <- subset(objdf, Time >= from)
  		}
  		if (!is.null(to)){
  			objdf <- subset(objdf, Time <= to)
  		}
	} else if (qi == "Hazard Ratio"){
	  	objdf <- data.frame(obj$Xj, obj$HR, obj$Comparison)
	  	names(objdf) <- c("Xj", "HR", "Comparison")
	} else if (qi == "Relative Hazard"){
	  	spalette <- NULL
	  	objdf <- data.frame(obj$Xj, obj$HR)
	  	names(objdf) <- c("Xj", "HR")
	} else if (qi == "First Difference"){
		spalette <- NULL
		objdf <- data.frame(obj$Xj, obj$FirstDiff, obj$Comparison)
		names(objdf) <- c("Xj", "FirstDiff", "Comparison")
	}

	# Plot
	  if (qi == "Hazard Rate"){
	  	if (!is.null(objdf$strata)) {
	  		ggplot(objdf, aes(x = Time, y = HRate, colour = factor(Comparison))) +
	  			geom_point(alpha = I(palpha), size = psize) +
	  			geom_smooth(method = smoother, size = lsize, se = 	FALSE) +
	  			facet_grid(.~ Strata) +
	  			scale_y_continuous()+
	  			scale_x_continuous() +
	  			xlab(xlab) + ylab(ylab) +
	  			scale_colour_brewer(palette = spalette, name = leg.name) +
	  			ggtitle(title) +
	  			guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	  		theme_bw(base_size = 15)
    	} else if (is.null(objdf$strata)){
	      	ggplot(objdf, aes(Time, HRate, colour = factor(Comparison))) +
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
	} else if (qi == "Relative Hazard"){
		ggplot(objdf, aes(Xj, HR)) +
		    geom_point(shape = 21, alpha = I(palpha), size = psize, colour = colour) +
		    geom_hline(aes(yintercept = 1), linetype = "dotted") +
		    scale_y_continuous() +
			scale_y_continuous() +    
		    xlab(xlab) + ylab(ylab) +
		    ggtitle(title) +
		    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		    theme_bw(base_size = 15)
	} else if (qi == "First Difference"){
    	ggplot(objdf, aes(Xj, FirstDiff, group = Comparison)) +
        	geom_point(shape = 21, alpha = I(palpha), size = psize, colour = colour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE) +
	        geom_hline(aes(yintercept = 0), linetype = "dotted") +
	        scale_y_continuous()+
	        scale_x_continuous() +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
	} else if (qi == "Hazard Ratio"){
		ggplot(objdf, aes(Xj, HR)) +
	        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = colour) +
	        stat_smooth(method = smoother, size = lsize, colour = colour, se = FALSE) +
	        geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        scale_y_continuous()+
	        scale_x_continuous() +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
    }
}