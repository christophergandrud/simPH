#' Simulate non-linear hazards for a range of values
#'
#' \code{coxsimPoly} simulates hazards for polynomial covariate effects.
#' @param obj a coxph fitted model object with a polynomial coefficient.
#' @param b character string name of the coefficient you would like to simulate.
#' @param pow numeric polynomial used in \code{coxph}.  
#' @param X numeric vector of values of X to simulate relative hazards for.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @description Simulates hazards for polynomial covariate effects.
#' 
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + 
#'              acutediz + hosp01  + hhosleng + mandiz01 + femdiz01 + 
#'              peddiz01 + orphdum + natreg + I(natreg^2) + vandavg3 + 
#'              wpnoavg3 + condavg3 + orderent + stafcder, 
#'             data = CarpenterFdaData)
#'
#' @references Keele, Luke. 2010. “Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.” Political Analysis 18(2): 189–205.
#'
#' Carpenter, Daniel P. 2002. “Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.” American Journal of Political Science 46(3): 490–505.
#' @import MSBVAR plyr reshape2 survival
#' @export 

coxsimPoly <- function(obj, b, pow = 2, X = 1, nsim = 1000) 
{
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	    
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	bpos <- match(b, dfn)

	NamesLoc <- function(p){
		Temp <- paste0("I.", b, ".", p, ".")
		match(Temp, dfn)
  	}
	pows <- 2:pow
  	NamePow <- sapply(pows, NamesLoc, simplify = TRUE)

	NamePow <- c(bpos, NamePow)
  	Drawn <- data.frame(Drawn[, NamePow])
  	Drawn$ID <- 1:nsim

}