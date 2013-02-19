#' Plot simulated penalised spline hazards.
#'
#' \code{ggspline} uses ggplot2 and scatter3d to plot quantities of interest from \code{simspline} objects, including relative hazards, first differences, hazard ratios, and hazard rates.
#'
#' @param obj a simlinear object
#' @param qi character string indicating what quantity of interest you would like to calculate. Can be \code{'Relative Hazard'}, \code{'First Difference'}, \code{'Hazard Ratio'}, or \code{'Hazard Rate'}. Default is \code{qi = 'Relative Hazard'}. 
#' @param FacetTime a numeric vector of points in time where you would like to plot Hazard Rates in a facet grid. Only relevant if \code{qi == 'Hazard Rate'}. Note: the values of Facet Time must exactly match values of the \code{time} element of \code{obj}.
#' @param xlab a label for the plot's x-axis.
#' @param ylab a label of the plot's y-axis. The default uses the value of \code{qi}.
#' @param zlab a label for the plot's z-axis. Only relevant if \code{qi = "Hazard Rate"} and \code{FacetTime == NULL}.
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
#' @param surface plot surface. Default is \code{surface = TRUE}. Only relevant if \code{qi == 'Relative Hazard'} and \code{FacetTime = Null}.
#' @param fit one or more of \code{"linear"}, \code{"quadratic"}, \code{"smooth"}, \code{"additive"}; to display fitted surface(s); partial matching is supported – e.g., \code{c("lin", "quad")}. Only relevant if \code{qi == 'Relative Hazard'} and \code{FacetTime = Null}.
#' @param ... other arguments passed to specific methods
#' @return a ggplot2 object. See \code{\link{scatter3d}} for values from \code{scatter3d} calls.
#'  
#' @description Uses \code{ggplot2} and \code{scatter3d} to plot the quantities of interest from \code{simspline} objects, including relative hazards, first differences, hazard ratios, and hazard rates. If currently does not support hazard rates for multiple strata.
#'
#' It can graph hazard rates as a 3D plot using \code{\link{scatter3d}} with the dimensions: Time, Hazard Rate, and the value of Xj. You can also choose to plot hazard rates for a range of values of Xj in two dimensional plots at specific points in time. Each plot is arranged in a facet grid.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#' 
#' # Load survival package
#' library(survival)
#' 
#' # Run basic model
#' # From Keele (2010) replication data
#' M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 + acutediz + hosp01  + pspline(hospdisc, df = 4) + pspline(hhosleng, df = 4) + mandiz01 + femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + wpnoavg3 + pspline(condavg3, df = 4) + pspline(orderent, df = 4) + pspline(stafcder, df = 4), data = CarpenterFdaData)
#'
#' # Simulate Relative Hazards for orderent
#' Sim1 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#'                        bdata = CarpenterFdaData$orderent,
#'                        qi = "Relative Hazard",
#'                        Xj = seq(1, 30, by = .1), ci = "90")
#'
#' # Simulate Hazard Rate for orderent
#' Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#'                     bdata = CarpenterFdaData$orderent,
#'                     qi = "Hazard Rate",
#'                     Xj = seq(1, 30, by = 2), ci = "90", nsim = 10)  
#'                        
#' # Plot relative hazard
#' ggspline(Sim1, xlab = "\n orderent", palpha = 1)
#' 
#' # 3D plot hazard rate
#' ggspline(Sim2, qi = "Hazard Rate", zlab = "orderent", fit = "quadratic")
#'
#' # Create a time grid plot
#' # Find all points in time where baseline hazard was found
#' unique(Sim2$time)
#' 
#' # Round time values so they can be exactly matched with FacetTime
#' Sim2$time <- round(Sim2$time, digits = 2)
#' 
#' # Create plot
#' ggspline(Sim2, qi = "Hazard Rate", FacetTime = c(6.21, 25.68, 100.64, 202.36))
#' 
#' @seealso \code{\link{coxsimLinear}}, \code{\link{ggtvc}},  \code{\link{ggplot2}}, and \code{\link{scatter3d}} 
#' 
#' 
#' 
#' 
#' @import ggplot2 car
#' @export

ggspline <- function(obj, qi = "Relative Hazard", FacetTime = NULL, from = NULL, to = NULL, xlab = NULL, ylab = NULL, zlab = NULL, title = NULL, smoother = "auto", spalette = "Set1", leg.name = "", lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3", psize = 1, palpha = 0.1, surface = TRUE, fit = "linear", ...)
{
	if (!inherits(obj, "simspline")){
    	stop("must be a simspline object")
    }
    # Create y-axis label
    if (is.null(ylab)){
    	ylab <- paste(qi, "\n")
    } else {
    	ylab <- ylab
    }
    # Subset simspline object & create a data frame of important variables
	if (qi == "Hazard Rate"){
		spallette <- NULL
		if (is.null(obj$strata)){
			objdf <- data.frame(obj$time, obj$HRate, obj$Xj)
			names(objdf) <- c("Time", "HRate", "Xj")
		} else if (!is.null(obj$strata)) {
		objdf <- data.frame(obj$time, obj$HRate, obj$strata, obj$Xj)
		names(objdf) <- c("Time", "HRate", "Strata", "Xj")
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
		objdf <- data.frame(obj$Xj, obj$HR, obj$Comparison)
		names(objdf) <- c("Xj", "FirstDiff", "Comparison")
	}

	# Plots
	if (qi == "Relative Hazard"){
		ggplot(objdf, aes(Xj, HR)) +
		    geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
		    geom_hline(aes(yintercept = 1), linetype = "dotted") +
		    scale_y_continuous() +
			scale_x_continuous() +    
		    xlab(xlab) + ylab(ylab) +
		    ggtitle(title) +
		    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
		    theme_bw(base_size = 15)
	} else if (qi == "First Difference"){
    	ggplot(objdf, aes(Xj, FirstDiff, group = Comparison)) +
        	geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 0), linetype = "dotted") +
	        scale_y_continuous()+
	        scale_x_continuous() +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
	} else if (qi == "Hazard Ratio"){
		ggplot(objdf, aes(Xj, HR)) +
	        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        scale_y_continuous()+
	        scale_x_continuous() +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
    } else if (qi == "Hazard Rate" & is.null(FacetTime)){
    	with(objdf, scatter3d(x = Time, y = HRate, z = Xj,
    						  xlab = xlab, ylab = ylab, zlab = zlab,
    						  surface = surface,
    						  fit = fit))
    } else if (qi == "Hazard Rate" & !is.null(FacetTime)){
		SubsetTime <- function(f){
		  CombObjdf <- data.frame()
		  for (i in f){
		    Temps <- objdf
		    TempsSub <- subset(Temps, Time == i)
		    CombObjdf <- rbind(CombObjdf, TempsSub)
		  }
		  CombObjdf
		}
		objdfSub <- SubsetTime(FacetTime)
		ggplot(objdfSub, aes(Xj, HRate)) +
	        geom_point(shape = 21, alpha = I(palpha), size = psize, colour = pcolour) +
	        geom_smooth(method = smoother, size = lsize, se = FALSE, color = lcolour) +
	        geom_hline(aes(yintercept = 1), linetype = "dotted") +
	        facet_grid(.~Time) +
	        scale_y_continuous()+
	        scale_x_continuous() +
	        xlab(xlab) + ylab(ylab) +
	        ggtitle(title) +
	        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	        theme_bw(base_size = 15)
    }
}