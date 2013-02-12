#' Simulate quantities of interest for linear multiplicative interactions.
#'
#' \code{coxsimInteract} simulates quantities of interest for linear multiplicative interactions.
#' @param obj a coxph fitted model object with a linear multiplicative interaction.
#' @param bs character vector of the two constitutive variables' names.
#' @param qi quantities of interest to simulate. Values can be \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. If \code{qi = "Hazard Rate"} and the \code{coxph} model has strata, then hazard rates for each strata will also be calculated.
#' @X1 numeric vector of fitted values of \code{X1} to simulate for.3
#' @X2 numeric vector of fitted values of \code{X2} to simulate for.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#'
#' @references Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. “Understanding Interaction Models: Improving Empirical Analyses.” Political Analysis 14(1): 63–82.
#'
#' @returm a siminteract class object
#' @import MSBVAR plyr reshape2 survival
#' @export

coxsimInteract <- function(obj, bs, qi = "Relative Hazard", X1 =NULL, X2 = NULL, nsim = 1000, ci = "95")
{
	# Parameter estimates & Variance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)
	
	# Draw covariate estimates from the multivariate normal distribution	    
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset data frame to only include interaction constitutive terms and 
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
			Simb$ME <- exp(Simb[, 1] + (Simb[, 3] * Simb[, 4]))	
		}
	}

	else if (qi == "Relative Hazard"){
		X2 <- rep(0, length(X1))
		Xs <- data.frame(X1, X2)
	Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	Simb <- merge(Simb, Xs)
		Simb$HR <- exp((Simb$X1 - Simb$X2) * Simb$Coef)	
	} 
	else if (qi == "First Difference"){
		if (length(X1) != length(X2)){
	  stop("X1 and X2 must be the same length.")
	} 
	else {
	    Xs <- data.frame(X1, X2)
	    Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
	    Simb$FirstDiff <- (exp((Simb$X1 - Simb$X2) * Simb$Coef) - 1) * 100
		}
	}
	else if (qi == "Hazard Ratio"){
	if (length(X2) > 1 & length(X1) != length(X2)){
	  stop("X1 and X2 must be the same length.")
	}
	else {
		if (length(X1) > 1 & length(X2) == 1){
			X2 <- rep(0, length(X1))
		}
		Xs <- data.frame(X1, X2)
		Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		    Simb$HR <- exp((Simb$X1 - Simb$X2) * Simb$Coef)	
	} 
	}
	else if (qi == "Hazard Rate"){
	if (length(X2) > 1 & length(X1) != length(X2)){
	  stop("X1 and X2 must be the same length.")
	}
	else {
		if (length(X1) > 1 & length(X2) == 1){
			X2 <- rep(0, length(X1))
		}
		Xs <- data.frame(X1, X2)
		Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
	    Simb <- merge(Simb, Xs)
		    Simb$HR <- exp((Simb$X1 - Simb$X2) * Simb$Coef)	 
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
		SubVar <- "X1"
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
	class(SimbPerc) <- "siminteract"
	SimbPerc
}