#' Constrict simulations to a defined interval.
#'
#' \code{IntervalConstrict} is an internal function to constrict a set of simulations to user defined interval.
#'
#' @param Simb character string naming the data frame with the simulations.
#' @param SubVar character vector the variable names to subset the simulations by.
#' @param QI character string labeling the quantitiy of interest.
#' @param spin logical for whether or not to use the shortest probability interval or the central interval.
#' @param ci numeric confidence interval measure.
#' 
#' @importFrom plyr ddply
#' @export

IntervalConstrict <- function(Simb = Simb, SubVar = SubVar, QI = QI, spin = FALSE, ci = 0.95)
{
	
}
