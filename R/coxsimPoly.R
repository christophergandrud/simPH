#' Simulate non-linear hazards for a range of values
#'
#' \code{coxsimPoly} simulates hazards for polynomial covariate effects.
#' @param obj a coxph fitted model object with a polynomial coefficient.
#' @param b character string name of the coefficient you would like to simulate.
#' @param pow numeric polynomial used in \code{coxph}.  
#' @param X numeric vector of values of X to simulate relative hazards for.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#' @return a simpoly class object.
#' @description Simulates relative hazards for polynomial covariate effects.
#'
#' Note, you must use \code{\link{I}} to create the polynomials.
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
#' # Simulate simpoly class object
#' # Sim1 <- coxsimPoly(M1, b = "natreg", pow = 3, X = seq(1, 150, by = 5))
#' 
#' @references Keele, Luke. 2010. “Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.” Political Analysis 18(2): 189–205.
#'
#' Carpenter, Daniel P. 2002. “Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.” American Journal of Political Science 46(3): 490–505.
#' @seealso \code{\link{ggpoly}}, \code{\link{rmultinorm}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#' @import MSBVAR plyr reshape2 survival
#' @export 

coxsimPoly <- function(obj, b, pow = 2, X, nsim = 1000, ci = "95") 
{
	# Parameter estimates & Variance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	  
	# Draw covariate estimates from the multivariate normal distribution	   
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset data frame to only include polynomial constitutive terms.
	bpos <- match(b, dfn)
	NamesLoc <- function(p){
		Temp <- paste0("I.", b, ".", p, ".")
		match(Temp, dfn)
  	}
	pows <- as.numeric(2:pow)
  	NamePow <- sapply(pows, NamesLoc, simplify = TRUE)

	NamePow <- c(bpos, NamePow)
  	Drawn <- data.frame(Drawn[, NamePow])
  	VNames <- names(Drawn)
	powFull <- as.numeric(1:pow)

	# Function to Multiply covariates by polynomials
  	Fitted <- function(VN, x, p){
  		Temp <- outer(Drawn[, VN], x^p)
  		TempDF <- data.frame(melt(Temp))
  		TempDF <- TempDF[, "value"]
  		TempDF
  	}

  	# Create combined relative hazards
    CombinedDF <- data.frame()
  	for (i in X){
  		TempComb <- mapply(Fitted, VN = VNames, x = i, p = powFull)
  		TempComb <- data.frame(TempComb)
  		TempComb$X <- i
  		TempComb$ID <- 1:nsim
  		CombinedDF <- rbind(CombinedDF, TempComb)
  	}
  	CombinedDF$RH <- exp(rowSums(CombinedDF[, VNames]))

  	# Drop simulations outside of 'confidence bounds'
	if (ci == "all"){
	    PolySimPerc <- CombinedDF 
	  } else if (ci == "95"){
	    PolySimPerc <- ddply(CombinedDF, .(X), transform, Lower = RH < quantile(RH, c(0.025)))
	    PolySimPerc <- ddply(PolySimPerc, .(X), transform, Upper = RH > quantile(RH, 0.975))
	  } else if (ci == "90"){
	    PolySimPerc <- ddply(CombinedDF, .(X), transform, Lower = RH < quantile(RH, c(0.05)))
	    PolySimPerc <- ddply(PolySimPerc, .(X), transform, Upper = RH > quantile(RH, 0.95))
	  } else if (ci == "99"){
	    PolySimPerc <- ddply(CombinedDF, .(X), transform, Lower = RH < quantile(RH, c(0.005)))
	    PolySimPerc <- ddply(PolySimPerc, .(X), transform, Upper = RH > quantile(RH, 0.995))
	  }
  	if (ci != "all"){
    	PolySimPerc <- subset(PolySimPerc, Lower == FALSE & Upper == FALSE)
  	}

	# Clean up
	PolySimPerc <- PolySimPerc[, c("ID", "X", "RH")]
	class(PolySimPerc) <- "simpoly"
	PolySimPerc
}





