#' Simulate quantities of interest for a range of values for a polynomial nonlinear effect from Cox Proportional Hazards models
#'
#' \code{coxsimPoly} simulates quantities of interest for polynomial covariate effects estimated from Cox Proportional Hazards models. These can be plotted with \code{\link{simGG}}.
#' @param obj a \code{\link{coxph}} class fitted model object with a polynomial coefficient. These can be plotted with \code{\link{simGG}}.
#' @param b character string name of the coefficient you would like to simulate.
#' @param qi quantity of interest to simulate. Values can be \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. If \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard rates for each strata will also be calculated.
#' @param pow numeric polynomial used in \code{coxph}.  
#' @param Xj numeric vector of fitted values for \code{b} to simulate for.
#' @param Xl numeric vector of values to compare \code{Xj} to. If \code{NULL}, then it is authomatically set to 0.
#' @param nsim the number of simulations to run per value of \code{Xj}. Default is \code{nsim = 1000}.
#' @param ci the proportion of simulations to keep. The default is \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE} then \code{ci} is the confidence level of the shortest probability interval. Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest probability interval rather than the middle simulations.
#'
#' @return a \code{simpoly} class object.
#' @details Simulates quantities of interest for polynomial covariate effects. For example if a nonlinear effect is modeled with a second order polynomial--i.e. \eqn{\beta_{1}x_{i} + \beta_{2}x_{i}^{2}}--we can once again draw \eqn{n} simulations from the multivariate normal distribution for both \eqn{\beta_{1}} and \eqn{\beta_{2}}. Then we simply calculate quantities of interest for a range of values and plot the results as before. For example, we find the first difference for a second order polynomial with:
#' \deqn{\%\triangle h_{i}(t) = (\mathrm{e}^{\beta_{1}x_{j-1} + \beta_{2}x_{j-l}^{2}} - 1) * 100}

#' where \eqn{x_{j-l} = x_{j} - x_{l}}.
#'
#' Note, you must use \code{\link{I}} to create the polynomials.
#' 
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz +
#' 				hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum + 
#' 				natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 + 
#'				condavg3 + orderent + stafcder, data = CarpenterFdaData)
#' 
#' # Simulate simpoly First Difference
#' Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference", 
#'						pow = 3, Xj = seq(1, 150, by = 5))
#'
#' # Simulate simpoly Hazard Ratio with spin probibility interval
#' Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio", 
#'						pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)
#' 
#' @references Keele, Luke. 2010. ''Proportionally Difficult: Testing for Nonproportional Hazards in Cox Models.'' Political Analysis 18(2): 189-205.
#'
#' Carpenter, Daniel P. 2002. ''Groups, the Media, Agency Waiting Costs, and FDA Drug Approval.'' American Journal of Political Science 46(3): 490-505.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' American Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient Shortest Probability Intervals.'' Arvix. \url{http://arxiv.org/pdf/1302.2142v1.pdf}.
#'
#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#' @importFrom reshape2 melt
#' @importFrom MSBVAR rmultnorm
#' @importFrom survival basehaz
#' @importFrom plyr rename
#' @export 

coxsimPoly <- function(obj, b, qi = "Relative Hazard", pow = 2, Xj = NULL, Xl = NULL, nsim = 1000, ci = 0.95, spin = FALSE) 
{
  QI <- NULL
	# Ensure that qi is valid
	qiOpts <- c("Relative Hazard", "First Difference", "Hazard Rate", "Hazard Ratio")
	TestqiOpts <- qi %in% qiOpts
	if (!isTRUE(TestqiOpts)){
  	stop("Invalid qi type. qi must be 'Relative Hazard', 'Hazard Rate', 'First Difference', or 'Hazard Ratio'")
	}

	# Find X_{jl}
	if (length(Xj) != length(Xl) & !is.null(Xl)){
		stop("Xj and Xl must be the same length.")
	}
	if (is.null(Xl) & qi != "Hazard Rate") {
		message("All Xl set at 0.")
		Xjl <- Xj
  } else if (!is.null(Xl) & qi == "Relative Hazard") {
    message("All Xl set to 0.")
    Xjl <- Xl
	} else {
  	Xbound <- cbind(Xj, Xl)
  	Xjl <- Xbound[, 1] - Xbound[, 2]
	}

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

    Simb <- data.frame()
  	for (i in Xjl){
  		TempComb <- mapply(Fitted, VN = VNames, x = i, p = powFull)
  		TempComb <- data.frame(TempComb)
  		TempComb$Xjl <- i
  		TempComb$ID <- 1:nsim
  		Simb <- rbind(Simb, TempComb)
  	}

   	# Create combined quantities of interest
    if (qi == "Relative Hazard"){
      Simb$QI <- exp(rowSums(Simb[, VNames]))
  	} else if (qi == "Hazard Ratio"){
  		Simb$QI <- exp(rowSums(Simb[, VNames]))
  	}
  	else if (qi == "First Difference"){
  		Simb$QI <- (exp(rowSums(Simb[, VNames])) - 1) * 100
  	}
  	else if (qi == "Hazard Rate"){
  		Simb$HR <- exp(rowSums(Simb[, VNames]))
  	  	bfit <- basehaz(obj)
  	  	bfit$FakeID <- 1
  	  	Simb$FakeID <- 1
        bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
        SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
        SimbCombDT <- SimbDT[bfitDT, allow.cartesian=TRUE]
        Simb <- data.frame(SimbCombDT)
  	  	Simb$QI <- Simb$hazard * Simb$HR 
  	  	Simb <- Simb[, -1]
  	}

  	# Drop simulations outside of 'confidence bounds'
	if (qi != "Hazard Rate"){
		SubVar <- "Xjl"
	} else if (qi == "Hazard Rate"){
		Simb <- rename(Simb, replace = c("Xjl" = "HRValue"))
		SubVar <- c("time", "HRValue")
	}

  SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar, qi = qi,
                                QI = QI, spin = spin, ci = ci)

	# Clean up
	class(SimbPerc) <- c("simpoly", qi)
	SimbPerc
}