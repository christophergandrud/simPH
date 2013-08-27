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
	Xj <- QI <- Time <- HRValue <- HRate <- Strata <- NULL
	if (isTRUE(hr) & !isTRUE(strata)){
		byVars <- c("Time", "HRValue")
	}
	else if (isTRUE(hr) & !isTRUE(strata)){
		byVars <- c("Time", "HRValue", "Strata")
	}
	if (!isTRUE(hr)){
		Linesdf <- ddply(df, byVars, transform, Median = median(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Max = max(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Min = min(QI))
		Linesdf <- ddply(Linesdf, byVars, transform, Lower50 = quantile(QI, 0.25))
		Linesdf <- ddply(Linesdf, byVars, transform, Upper50 = quantile(QI, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, byVars]), ] 
	}
	else if (isTRUE(hr) & !isTRUE(strata)){
		Linesdf <- ddply(df, byVars, transform, Median = median(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Max = max(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Min = min(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Lower50 = quantile(HRate, 0.25))
		Linesdf <- ddply(Linesdf, byVars, transform, Upper50 = quantile(HRate, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, c(1, 3)]), ] 

	}
	else if (isTRUE(hr) & isTRUE(strata)){
		Linesdf <- ddply(df, byVars, transform, Median = median(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Max = max(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Min = min(HRate))
		Linesdf <- ddply(Linesdf, byVars, transform, Lower50 = quantile(HRate, 0.25))
		Linesdf <- ddply(Linesdf, byVars, transform, Upper50 = quantile(HRate, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[, c("Time", "HRValue", "Strata")]), ] 		
	}
	return(Linesdf)
}