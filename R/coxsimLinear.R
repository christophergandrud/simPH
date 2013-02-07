#' Simulate hazards for linear time-constant covariates from Cox Proportional Hazards models.
#'
#' \code{coxsimLinear} simulates relative hazards, first differences, and hazard ratios for time-constant covariates from models estimated with \code{\link{coxph}} using the multivariate normal distribution.
#' @param obj a coxph fitted model object.
#' @param b character string name of the coefficient you would like to simulate.
#' @param qi quantity of interest to simulate. Values can be \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{"qi = "Relative Hazard"}. If code{"qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard rates for each strata will also be calculated.
#' @param Xj numeric vector of values of X to simulate for.
#' @param Xl numeric vector of values to compare \code{Xj} to. Note if \code{qi = "Relative Hazard"} only \code{Xj} is relevant.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#' @return a simlinear object
#'
#'
#' @description Simulates relative hazards, first differences, hazard ratios, and hazard rates for linear time-constant covariates from Cox Proportional Hazard models.
#'
#' @examples
#' # Load Carpenter (2002) data
#' # data("CarpenterFdaData")
#'
#' # Load survival package
#' # library(survival)
#'
#' # Run basic model
#' # M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum +natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + condavg3 + orderent + stafcder, data = CarpenterFdaData)
#'
#' # Simulate values
#' # Sim1 <- coxsimLinear(M1, b = "stafcder", qi = "Hazard Ratio", Xj = c(1237, 1600), Xl = c(1000, 1000), ci = "99")
#'
#' @seealso \code{\link{ggLinear}}, \code{\link{rmultinorm}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#' @references Licht, Amanda A. 2011. “Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.” Political Analysis 19: 227–43.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. “Making the Most of Statistical Analyses: Improving Interpretation and Presentation.” American Journal of Political Science 44(2): 347–61.
#' @import MSBVAR plyr reshape2 survival
#' @export

coxsimLinear <- function(obj, b, qi = "Relative Hazard", Xj = 1, Xl = 0, nsim = 1000, ci = "95")
{	
	# Parameter estimates & Varance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	
	# Draw covariate estimates from the multivariate normal distribution	    
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset simulations to only include b
	bpos <- match(b, dfn)
	Simb <- data.frame(DrawnDF[, bpos])
	names(Simb) <- "Coef"

  # Find quantity of interest
  if (qi == "Relative Hazard"){
  	Xl <- rep(0, length(Xj))
  	Xs <- data.frame(Xj, Xl)
    Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	Simb <- merge(Simb, Xs)
  	Simb$HR <- exp((Simb$Xj - Simb$Xl) * Simb$Coef)	
  } 
  else if (qi == "First Difference"){
  	if (length(Xj) != length(Xl)){
      stop("Xj and Xl must be the same length.")
    } 
    else {
	    Xs <- data.frame(Xj, Xl)
	    Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
	    Simb$FirstDiff <- (exp((Simb$Xj - Simb$Xl) * Simb$Coef) - 1) * 100
  	}
  }
  else if (qi == "Hazard Ratio"){
    if (length(Xl) > 1 & length(Xj) != length(Xl)){
      stop("Xj and Xl must be the same length.")
    }
    else {
    	if (length(Xj) > 1 & length(Xl) == 1){
    		Xl <- rep(0, length(Xj))
    	}
    	Xs <- data.frame(Xj, Xl)
    	Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
  	    Simb$HR <- exp((Simb$Xj - Simb$Xl) * Simb$Coef)	
    } 
  }
  else if (qi == "Hazard Rate"){
    if (length(Xl) > 1 & length(Xj) != length(Xl)){
      stop("Xj and Xl must be the same length.")
    }
    else {
    	if (length(Xj) > 1 & length(Xl) == 1){
    		Xl <- rep(0, length(Xj))
    	}
    	Xs <- data.frame(Xj, Xl)
    	Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
  	    Simb$HR <- exp((Simb$Xj - Simb$Xl) * Simb$Coef)	 
	  	bfit <- basehaz(obj)
	  	bfit$FakeID <- 1
	  	Simb$FakeID <- 1
	  	Simb <- merge(bfit, Simb, by = "FakeID")
	  	Simb$HRate <- Simb$hazard * Simb$HR 
	  	Simb <- Simb[, -1]
  	}
  }

  # Drop simulations outside of 'confidence bounds'
  if (qi != "Hazard Rate"){
  	SubVar <- "Xj"
  } else if (qi == "Hazard Rate"){
  	SubVar <- "time"
  }
  if (ci == "all"){
    SimbPerc <- Simb 
  } else if (ci == "90"){
    SimbPerc <- ddply(Simb, SubVar, transform, Lower = HR < quantile(HR, c(0.05)))
    SimbPerc <- ddply(SimbPerc, SubVar, transform, Upper = HR > quantile(HR, 0.95))
  } else if (ci == "95"){
    SimbPerc <- ddply(Simb, SubVar, transform, Lower = HR < quantile(HR, c(0.025)))
    SimbPerc <- ddply(SimbPerc, SubVar, transform, Upper = HR > quantile(HR, 0.975))
  } else if (ci == "99"){
    SimbPerc <- ddply(Simb, SubVar, transform, Lower = HR < quantile(HR, c(0.005)))
    SimbPerc <- ddply(SimbPerc, SubVar, transform, Upper = HR > quantile(HR, 0.995))
  }
  if (ci != "all"){
    SimbPerc <- subset(SimbPerc, Lower == FALSE & Upper == FALSE)
  }  

  # Final clean up
  class(SimbPerc) <- "simlinear"
  SimbPerc

}