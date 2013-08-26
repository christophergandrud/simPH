#' Transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent
#'
#' \code{MinMaxLines} is an internal function to transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent.
#'
#' @param df the data frame created by the simGG method from a simulation class object.
#' @param hr logical indicating whether or not \code{df} contains a hazard rate.
#' @param strata logical indicating whether or not \code{df} contains a stratified hazard rate. 
#'
#' @importFrom plyr ddply mutate
#' @keywords internals
#' @noRd

MinMaxLines <- function(df, hr = FALSE, strata = FALSE){
	if (!isTRUE(hr)){
		Linesdf <- ddply(df, .(Xj), transform, Median = median(QI))
		Linesdf <- ddply(Linesdf, .(Xj), transform, Max = max(QI))
		Linesdf <- ddply(Linesdf, .(Xj), transform, Min = min(QI))
		Linesdf <- ddply(Linesdf, .(Xj), transform, Lower50 = quantile(QI, 0.25))
		Linesdf <- ddply(Linesdf, .(Xj), transform, Upper50 = quantile(QI, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[[1]]), ] 
		Linesdf <- Linesdf[, -2]
	}
	else if (isTRUE(isTRUE(hr) & !isTRUE(strata))){
		Linesdf <- ddply(df, .(Time), transform, Median = median(HRate))
		Linesdf <- ddply(Linesdf, .(Time), transform, Max = max(HRate))
		Linesdf <- ddply(Linesdf, .(Time), transform, Min = min(HRate))
		Linesdf <- ddply(Linesdf, .(Time), transform, Lower50 = quantile(HRate, 0.25))
		Linesdf <- ddply(Linesdf, .(Time), transform, Upper50 = quantile(HRate, 0.75))

		Linesdf <- Linesdf[!duplicated(Linesdf[[1]]), ] 
		# Linesdf <- Linesdf[, -2]		
	}
}