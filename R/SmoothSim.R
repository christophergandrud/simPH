#' Smooth values from one simulation using smoothing splines
#'
#' @param x a vector containing the x value.
#' @param y a vector containing the y value.
#' @param df numeric. The desired equivalent number of degrees of freedom (trace of the smoother matrix).
#'
#' @references Chambers, J. M. and Hastie, T. J. (1992) Statistical Models in S, Wadsworth & Brooks/Cole.
#'
#' Green, P. J. and Silverman, B. W. (1994) Nonparametric Regression and Generalized Linear Models: A Roughness Penalty Approach. Chapman and Hall.
#'
#'Hastie, T. J. and Tibshirani, R. J. (1990) Generalized Additive Models. Chapman and Hall.
#'
#' @importFrom stats smooth.spline
#' @keywords internals
#' @noRd

SmoothOneSim <- function(x, y, df = 10){
  TempXY <- cbind(x, y)
  TempOut <- smooth.spline(TempXY, df = df)
  Out <- fitted(TempOut)
  return(Out)
}

#' Smooth values for all simulations
#'
#' @param SimIn data frame. Pre-smoothed simulation
#' @param xaxis character string. The column that will form the xaxis in the plot.  
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#'
#' @keywords internals
#' @noRd
SmoothSimulations <- function(SimIn, xaxis = "Xj"){
	# Drop simulations that do not have all values within
  # the central interval
	SimIn$dummy <- 1
	Sims <- dplyr::group_by(SimIn, SimID) 
 	Sims <- dplyr::mutate(Sims, Rows = sum(dummy))
  MaxRows <- max(Sims$Rows)
  Sims <- as.data.frame(Sims)
  Sims <- subset(Sims, Rows == MaxRows)
	names(Sims)[names(Sims) == xaxis] <- 'Xj'

  # Spline smooth
  Sims <- dplyr::group_by(Sims, SimID) 
  SimsFitted <- dplyr::mutate(Sims, QI = SmoothOneSim(Xj, QI))
	names(Sims)[names(Sims) == 'Xj'] <- xaxis
  SimsFitted
}