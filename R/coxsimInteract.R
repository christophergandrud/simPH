#' Simulate quantities of interest for linear multiplicative interactions.
#'
#' \code{coxsimInteract} simulates quantities of interest for linear multiplicative interactions.
#' @param obj a coxph fitted model object with a linear multiplicative interaction.
#' @param b1 character string of the first constitutive variable's name. Note \code{b1} and \code{b2} must be entered in the order in which they are entered into the \code{coxph} model.
#' @param b2 character string of the second constitutive variable's name.
#' @param qi quantities of interest to simulate. Values can be \code{"Marginal Effect"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. If \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard rates for each strata will also be calculated.
#' @param X1 numeric vector of fitted values of \code{b1} to simulate for. If \code{qi = "Marginal Effect"} then only \code{X2} can be set.
#' @param X2 numeric vector of fitted values of \code{b2} to simulate for. 
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#' @description Simulates marginal effects, first differences, hazard ratios, and hazard rates for linear multiplicative interactions. 
#'
#' Marginal effects are calculated as in Brambor et al. (2006) with the addition that we take the exponent, so that it resembles a hazard ratio. For an interaction between variables \eqn{X} and \eqn{Z} then the marginal effect for \eqn{X} is:
#' \deqn{ME_{X} = exp(\beta_{X} + \beta_{XZ}Z)}
#'
#' Unlike in \code{\link{coxsimtvc}} and \code{\link{coxsimLinear}} Hazard Ratios and First Differences can only be calculated by comparing a specified value of \code{X1} and \code{X2} to 0.
#'
#'
#' @examples 
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx, data = CarpenterFdaData)
#' 
#' # Simulate Marginal Effect of lethal for multiple values of prevgenx
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx", X2 = seq(2, 115, by = 2))
#'
#' @references Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. “Understanding Interaction Models: Improving Empirical Analyses.” Political Analysis 14(1): 63–82.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. “Making the Most of Statistical Analyses: Improving Interpretation and Presentation.” American Journal of Political Science 44(2): 347–61.
#'
#' @seealso \code{\link{gginteract}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}},
#' @return a siminteract class object
#' @import MSBVAR plyr reshape2 survival
#' @export

coxsimInteract <- function(obj, b1, b2, qi = "Marginal Effect", X1 = NULL, X2 = NULL, nsim = 1000, ci = "95")
{
	# Parameter estimates & Variance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	
	# Draw covariate estimates from the multivariate normal distribution	    
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset data frame to only include interaction constitutive terms and
	bs <- c(b1, b2) 
	bpos <- match(bs, dfn)
	binter <- paste0(bs[[1]], ".", bs[[2]])
	binter <- match(binter, dfn)
	NamesInt <- c(bpos, binter)
	Simb <- data.frame(Drawn[, NamesInt])

  # Find quantity of interest
	if (qi == "Marginal Effect"){
		if (!is.null(X1)){
			stop("For Marginal Effects only X2 should be specified.")
		} else{
			X2df <- data.frame(X2)
			names(X2df) <- c("X2")
			Simb <- merge(Simb, X2df)
			Simb$HR <- exp(Simb[, 1] + (Simb[, 3] * Simb[, 4])) 
		}
	}
	else if (qi == "First Difference"){
	  if (is.null(X1) | is.null(X2)){
	    stop("For First Differences both X1 and X2 should be specified.")
	  } else{
		Xs <- merge(X1, X2)
		names(Xs) <- c("X1", "X2")
		Xs$Comparison <- paste0(Xs[, 1], ", ", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		Simb$HR <- (exp((Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) + (Simb$X1 * Simb$X2 * Simb[, 3]) - 1) * 100)	
	  }
	}
	else if (qi == "Hazard Ratio"){
	  if (is.null(X1) | is.null(X2)){
	    stop("For Hazard Ratios both X1 and X2 should be specified.")
	  } else{
		Xs <- merge(X1, X2)
		names(Xs) <- c("X1", "X2")
		Xs$Comparison <- paste0(Xs[, 1], ", ", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		Simb$HR <- (exp((Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) + (Simb$X1 * Simb$X2 * Simb[, 3]) - 1) * 100)
	  }
	}
	else if (qi == "Hazard Rate"){
	  if (is.null(X1) | is.null(X2)){
	    stop("For Hazard Rates both X1 and X2 should be specified.")
	  } else{
		Xs <- merge(X1, X2)
		names(Xs) <- c("X1", "X2")
		Xs$HRValue <- paste0(Xs[, 1], ", ", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		Simb$HR <- exp((Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) + (Simb$X1 * Simb$X2 * Simb[, 3]))	 
	  	bfit <- basehaz(obj)
	  	bfit$FakeID <- 1
	  	Simb$FakeID <- 1
	  	Simb <- merge(bfit, Simb, by = "FakeID")
	  	Simb$HRate <- Simb$hazard * Simb$HR 
	  	Simb <- Simb[, -1]
	  }
	}

	# Drop simulations outside of 'confidence bounds'
	if (qi == "Relative Hazard" | qi == "First Difference" | qi == "Hazard Ratio"){
		SubVar <- "X1"
	} else if (qi == "Marginal Effect"){
		SubVar <- "X2"
	}else if (qi == "Hazard Rate"){
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
	class(SimbPerc) <- "siminteract"
	SimbPerc
}