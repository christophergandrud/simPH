#' Simulate hazards for linear time-constant covariates from Cox Proportional Hazards models.
#'
#' \code{coxsimLinear} simulates relative hazards, first differences, and hazard ratios for time-constant covariates from models estimated with \code{\link{coxph}}.
#' 
#' \code{obj} obj a coxph fitted model object.
#' @param b character string name of the coefficient you would like to simulate.
#' @param X numeric vector of values of X to simulate relative hazards for.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#'
#'
#'
#' @description Simulates relative hazards, first differences, hazard ratios, and hazard rates for linear time-constant covariates from Cox Proportional Hazard models.
#'
#'
#'
#' @import MSBVAR plyr reshape2 survival
#' @export

coxsimLinear <- function(obj, b, X, nsim = 1000, ci = "95")
{
	# Parameter estimates & Varance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	
	# Draw covaritate estiamtes from the multivariate normal distribution	    
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset simulations to only include b
	bpos <- match(b, dfn)

	# Function to Multiply covaraites by polynomials
  	Fitted <- function(bpos, x){
  		Temp <- outer(Drawn[, bpos], x)
  		TempDF <- data.frame(melt(Temp))
  		TempDF <- TempDF[, "value"]
  		TempDF
  	}

  	# Create combined relative hazards
    CombinedDF <- data.frame()
  	for (i in X){
  		TempComb <- mapply(Fitted, VN = VNames, x = i)
  		TempComb <- data.frame(TempComb)
  		TempComb$X <- i
  		TempComb$ID <- 1:nsim
  		CombinedDF <- rbind(CombinedDF, TempComb)
  	}
  	CombinedDF$RH <- exp(CombinedDF[, ])


}