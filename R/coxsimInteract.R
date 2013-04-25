#' Simulate quantities of interest for linear multiplicative interactions from \code{coxph} models.
#'
#' \code{coxsimInteract} simulates quantities of interest for linear multiplicative interactions using multivariate normal distributions.
#' @param obj a \code{coxph} fitted model object with a linear multiplicative interaction.
#' @param b1 character string of the first constitutive variable's name. Note \code{b1} and \code{b2} must be entered in the order in which they are entered into the \code{coxph} model.
#' @param b2 character string of the second constitutive variable's name.
#' @param qi quantities of interest to simulate. Values can be \code{"Marginal Effect"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. If \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard rates for each strata will also be calculated.
#' @param X1 numeric vector of fitted values of \code{b1} to simulate for. If \code{qi = "Marginal Effect"} then only \code{X2} can be set. If you want to plot the results, \code{X1} should have more than one value.
#' @param X2 numeric vector of fitted values of \code{b2} to simulate for. 
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE} then \code{ci} is the convidence level of the shortest probability interval. Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest proability interval rather than the middle simulations.
#'
#' @details Simulates marginal effects, first differences, hazard ratios, and hazard rates for linear multiplicative interactions. 
#' Marginal effects are calculated as in Brambor et al. (2006) with the addition that we take the exponent, so that it resembles a hazard ratio. For an interaction between variables \eqn{X} and \eqn{Z} then the marginal effect for \eqn{X} is:
#' \deqn{ME_{X} = exp(\beta_{X} + \beta_{XZ}Z)}
#'
#' Unlike in \code{\link{coxsimtvc}} and \code{\link{coxsimLinear}} Hazard Ratios and First Differences can only be calculated by comparing a specified value of \code{X1} and \code{X2} to 0.
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
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx", X2 = seq(2, 115, by = 2), spin = TRUE)
#'
#' @references Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. ''Understanding Interaction Models: Improving Empirical Analyses.'' Political Analysis 14(1): 63–82.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' American Journal of Political Science 44(2): 347–61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient Shortest Probablility Intervals.'' Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.
#'
#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}},
#' @return a siminteract class object
#' @import data.table
#' @importFrom plyr ddply
#' @importFrom survival basehaz
#' @importFrom MSBVAR rmultnorm
#' @export

coxsimInteract <- function(obj, b1, b2, qi = "Marginal Effect", X1 = NULL, X2 = NULL, nsim = 1000, ci = 0.95, spin = FALSE)
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
	  } else {
		Xs <- merge(X1, X2)
		names(Xs) <- c("X1", "X2")
		Xs$Comparison <- paste0(Xs[, 1], ", ", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		Simb$HR <- (exp((Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) + (Simb$X1 * Simb$X2 * Simb[, 3]) - 1) * 100)
	  }
	}
	else if (qi == "Hazard Rate"){
		Xl <- NULL
      message("Xl is ignored. All variables values other than b fitted at 0.") 
		Xs <- data.frame(Xj)
		Xs$HRValue <- paste(Xs[, 1])
	    Simb <- merge(Simb, Xs)
		Simb$HR <- exp((Simb$X1 * Simb[, 1]) + (Simb$X2 * Simb[, 2]) + (Simb$X1 * Simb$X2 * Simb[, 3]))	 
	  	bfit <- basehaz(obj)
	  	bfit$FakeID <- 1
	  	Simb$FakeID <- 1
		bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
		SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
		SimbCombDT <- SimbDT[bfitDT, allow.cartesian=TRUE]
	  	Simb$HRate <- Simb$hazard * Simb$HR 
	  	Simb <- Simb[, -1]
	}

	# Drop simulations outside of 'confidence bounds'
	if (qi == "Relative Hazard" | qi == "First Difference" | qi == "Hazard Ratio"){
		SubVar <- "X1"
	} else if (qi == "Marginal Effect"){
		SubVar <- "X2"
	}else if (qi == "Hazard Rate"){
		SubVar <- "time"
	}
	if (!isTRUE(spin)){
	    Bottom <- (1 - ci)/2
	    Top <- 1 - Bottom
	    SimbPerc <- eval(parse(text = paste0("ddply(Simb, SubVar, mutate, Lower = HR < quantile(HR,", 
	      Bottom, 
	      "))"
	    )))
	    SimbPerc <- eval(parse(text = paste0("ddply(SimbPerc, SubVar, mutate, Upper = HR > quantile(HR,", 
	      Top, 
	      "))"
	    )))
	}

	# Drop simulations outside of the shortest probability interval
	else if (isTRUE(spin)){
	    SimbPerc <- eval(parse(text = paste0("ddply(Simb, SubVar, mutate, Lower = HR < SpinBounds(HR, conf = ", ci, ", LowUp = 1))" )))
	    SimbPerc <- eval(parse(text = paste0("ddply(SimbPerc, SubVar, mutate, Upper = HR > SpinBounds(HR, conf = ", ci, ", LowUp = 2))" )))
	}

	SimbPerc <- subset(SimbPerc, Lower == FALSE & Upper == FALSE)

	# Final clean up
	class(SimbPerc) <- c("siminteract", qi)
	SimbPerc
}