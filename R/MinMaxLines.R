#' Transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent
#'
#' \code{MinMaxLines} is an internal function to transform the simulation object to include only the min, max, and the lower and upper bounds of the middle 50 percent.
#'
#' @param df the data frame created by the simGG method from a simulation class object
#' 
#' @importFrom plyr ddply mutate
#' @keywords internals
#' @noRd

MinMaxLines <- function(df){
	Linesdf <- ddply(df, .(Xj), transform, Median = median(QI))
	Linesdf <- ddply(Linesdf, .(Xj), transform, Max = max(QI))
	Linesdf <- ddply(Linesdf, .(Xj), transform, Min = min(QI))
	Linesdf <- ddply(Linesdf, .(Xj), transform, Lower50 = quantile(QI, 0.25))
	Linesdf <- ddply(Linesdf, .(Xj), transform, Upper50 = quantile(QI, 0.75))

	Linesdf <- Linesdf[!duplicated(Linesdf[[1]]), ] 
	Linesdf <- Linesdf[, -2]
}