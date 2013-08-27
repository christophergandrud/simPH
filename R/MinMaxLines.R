#' Transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent
#'
#' \code{MinMaxLines} is an internal function to transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent.
#'
#' @param df the data frame created by the simGG method from a simulation class object.
#' @param byVars character vector of the variables to subset the data frame by.
#' @param hr logical indicating whether or not \code{df} contains a hazard rate.
#' @param strata logical indicating whether or not \code{df} contains a stratified hazard rate. 
#'
#' @importFrom plyr ddply
#' @keywords internals
#' @noRd

MinMaxLines <- function(df, byVars = "Xj", hr = FALSE, strata = FALSE){
	Xj <- QI <- Time <- HRValue <- HRate <- Strata <- . <- NULL
	if (!isTRUE(hr)){
		Linesdf <- ddply(df, byVars, transform, Median = median(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Max = max(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Min = min(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Lower50 = quantile(QI, 0.25))
		Linesdf <- ddply(Linesdf, byVars, transform, Upper50 = quantile(QI, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, byVars]), ] 
		# Linesdf <- Linesdf[, -2]
	}
	else if (isTRUE(hr) & !isTRUE(strata)){
		Linesdf <- ddply(df, .(Time, HRValue), transform, Median = median(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue), transform, Max = max(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue), transform, Min = min(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue), transform, Lower50 = quantile(HRate, 0.25))
		Linesdf <- ddply(Linesdf, .(Time, HRValue), transform, Upper50 = quantile(HRate, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, c(1, 3)]), ] 

	}
	else if (isTRUE(hr) & isTRUE(strata)){
		Linesdf <- ddply(df, .(Time, HRValue, Strata), transform, Median = median(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue, Strata), transform, Max = max(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue, Strata), transform, Min = min(HRate))
		Linesdf <- ddply(Linesdf, .(Time, HRValue, Strata), transform, Lower50 = quantile(HRate, 0.25))
		Linesdf <- ddply(Linesdf, .(Time, HRValue, Strata), transform, Upper50 = quantile(HRate, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, c("Time", "HRValue", "Strata")]), ] 		
	}
	return(Linesdf)
}