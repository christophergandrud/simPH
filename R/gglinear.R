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

gglinear <- function(obj, qi = "Relative Hazard", from = NULL, to = NULL, xlab = NULL, ylab = NULL, title = NULL, xbreaks = NULL, xlabels = NULL, smoother = "auto", colour = "#A6CEE3", spalette = "Set1", leg.name = "", lsize = 2, psize = 1, palpha = 0.1, ...)
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
			objdf <- data.frame(obj$RealTime, obj$HRate, obj$Comparison)
			names(objdf) <- c("Time", "HRate", "Comparison")
		} else {
		objdf <- data.frame(obj$RealTime, obj$HRate, obj$strata, obj$Comparison)
		names(objdf) <- c("Time", "HRate", "Strata", "Comparison")
		}
		if (!is.null(from)){
			objdf <- subset(objdf, Time >= from)
  		}
  		if (!is.null(to)){
  			objdf <- subset(objdf, Time <= to)
  		}
	} else if (qi == "Hazard Ratio"){
	  	objdf <- data.frame(obj$RealTime, obj$HR, obj$Comparison)
	  	names(objdf) <- c("Time", "HR", "Comparison")
	} else if (qi == "Relative Hazard"){
	  	spalette <- NULL
	  	objdf <- data.frame(obj$RealTime, obj$HR)
	  	names(objdf) <- c("Time", "HR")
	} else if (qi == "First Difference"){
		spalette <- NULL
		objdf <- data.frame(obj$RealTime, obj$FirstDiff, obj$Comparison)
		names(objdf) <- c("Time", "FirstDiff", "Comparison")
	}
}