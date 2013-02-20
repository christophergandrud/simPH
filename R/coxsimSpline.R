#' Simulated quantities of interest for penalised splines from \code{coxph} models.
#'
#' \code{coxsimSpline} simulates quantities of interest from penalised splines using multivariate normal distributions.
#' @param obj a \code{coxph} fitted model object with a penalised spline.
#' @param bspline a character string of the full \code{\link{psline}} call used in \code{obj}.
#' @param bdata a numeric vector of splined variable's values.
#' @param qi quantity of interest to simulate. Values can be \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. Think carefully before using \code{qi = "Hazard Rate"}. You may be creating very many simulated values which can be very computationally intensive to do. Adjust the number of simulations per fitted value with \code{nsim}.
#' @param Xj numeric vector of values of X to simulate for.
#' @param Xl numeric vector of values to compare \code{Xj} to. Note if \code{qi = "Relative Hazard"} or \code{code = "Hazard"} only \code{Xj} is relevant.
#' @param nsim the number of simulations to run per value of X. Default is \code{nsim = 1000}.
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#' @return a simspline object
#'
#' @description Simulates relative hazards, first differences, hazard ratios, and hazard rates for penalised splines from Cox Proportional Hazards models. These can be plotted with \code{\link{ggspline}}. 
#'
#' Currently does not support simulating hazard rates form multiple stratified models.
#'
#' @examples
#' # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#' 
#' # Load survival package
#' library(survival)
#' 
#' # Run basic model
#' # From Keele (2010) replication data
#' M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 + acutediz + hosp01  + pspline(hospdisc, df = 4) + pspline(hhosleng, df = 4) + mandiz01 + femdiz01 + peddiz01 + orphdum + natreg + vandavg3 + wpnoavg3 + pspline(condavg3, df = 4) + pspline(orderent, df = 4) + pspline(stafcder, df = 4), data = CarpenterFdaData)
#'
#' # Simulate Relative Hazards for orderent
#' Sim1 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#'                        bdata = CarpenterFdaData$orderent,
#'                        qi = "Relative Hazard",
#'                        Xj = seq(2, 30, by = 3))
#'  
#' # Simulate Hazard Rates for orderent
#' Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#'                        bdata = CarpenterFdaData$orderent,
#'                        qi = "Hazard Rate",
#'                        Xj = seq(2, 53, by = 3),
#'                        nsim = 100)
#'
#'
#' @seealso \code{\link{gg}}, \code{\link{pspline}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#'
#' @references Luke Keele, "Replication data for: Proportionally Difficult: Testing for Nonproportional Hazards In Cox Models", 2010, http://hdl.handle.net/1902.1/17068 V1 [Version] 
#' 
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. “Making the Most of Statistical Analyses: Improving Interpretation and Presentation.” American Journal of Political Science 44(2): 347–61.
#' 
#' @import MSBVAR stringr reshape2 data.table
#' @export

coxsimSpline <- function(obj, bspline, bdata, qi = "Relative Hazard", Xj = 1, Xl = 0, nsim = 1000, ci = "95")
{ 
	if (nsim > 10 & qi == "Hazard Rate"){
		print(paste0("Warning: finding Hazard Rates with ", nsim, " simulations may take awhile.  Consider changing the number of simulations with nsim."))
	}
	# Find term number
	TermNum <- names(obj$pterms)
	bterm <- match(bspline, TermNum)
	if (is.na(bterm)){
		stop(paste0("Unable to find ", bspline, "."))
	}

	# Extract boundary knots for default Boundary.knots = range(x) & number of knots
	OA <- obj$assign
	ListKnots <- OA[bterm]
	NumKnots <- length(unlist(ListKnots))
	KnotIntervals <- levels(cut(bdata, breaks = NumKnots))

	# Parameter estimates & Variance/Covariance matrix
	Coef <- matrix(obj$coefficients)
	VC <- vcov(obj)

	# Draw covariate estimates from the multivariate normal distribution	   
	Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
	DrawnDF <- data.frame(Drawn)
	dfn <- names(DrawnDF)

	# Subset data frame to only spline variable coefficients.
	bword <-word(bspline, 1)
	b <- str_replace(bword, "pspline\\(", "")
	b <- str_replace(b, ",", "")

	NamesLoc <- function(p){
	  Temp <- paste0("ps.", b, ".", p)
	  match(Temp, dfn)
	}

	UpLim <- 2 + NumKnots
	CoeNum <- as.numeric(3:UpLim)
	NameCoe <- sapply(CoeNum, NamesLoc, simplify = TRUE)
	DrawnDF <- data.frame(DrawnDF[, NameCoe])

	# Match coefficients to knot interval
	IntervalStartAbs <- "\\(-?[0-9]*.[0-9]*e?\\+?[0-9]*,"
	IntervalFinishAbs <- ",-?[0-9]*.[0-9]*e?\\+?[0-9]*\\]" 
	IntervalStart <- str_match(KnotIntervals, IntervalStartAbs)
	IntervalStart <- str_replace(IntervalStart, "\\(", "")
	IntervalStart <- str_replace(IntervalStart, ",", "")
	IntervalStart <- as.numeric(IntervalStart)

	IntervalFinish <- str_match(KnotIntervals, IntervalFinishAbs)
	IntervalFinish <- str_replace(IntervalFinish, "\\]", "")
	IntervalFinish <- str_replace(IntervalFinish, ",", "")
	IntervalFinish <- as.numeric(IntervalFinish)

	CoefIntervals <- data.frame(names(DrawnDF), IntervalStart, IntervalFinish)
	names(CoefIntervals) <- c("CoefName", "IntervalStart", "IntervalFinish")

	# Melt Drawn DF to long format
	TempDF <- suppressMessages(data.frame(melt(DrawnDF)))
	names(TempDF) <- c("CoefName", "Coef")

	# Merge with CoefIntervals
	CoefIntervalsDT <- data.table(CoefIntervals, key = "CoefName")
	TempDT <- data.table(TempDF, key = "CoefName")
	TempCombDT <- TempDT[CoefIntervalsDT]
	TempDF <- data.frame(TempCombDT)

	# Merge in fitted X values
	MergeX <- function(f){
		CombinedDF <- data.frame()
		for (i in f){
		  Temps <- TempDF
		  Temps$X <- ifelse(TempDF[, 3] < i & i <= TempDF[, 4], i, NA)
		  Temps <- subset(Temps, !is.na(X))
		  CombinedDF <- rbind(CombinedDF, Temps)
		}
		CombinedDF
	}

	# Find quantities of interest
	if (qi == "Relative Hazard"){
		Xl <- NULL
		print("Xl is ignored")  
		Simb <- MergeX(Xj)
	    names(Simb) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
	    Simb$HR <- exp(Simb$Xj * Simb$Coef)	
	}
	else if (qi == "First Difference"){
	  	if (length(Xj) != length(Xl)){
	      stop("Xj and Xl must be the same length.")
	    } 
	    else {
		 	Xs <- data.frame(Xj, Xl)   	
			CombinedDF <- MergeX(Xj)
		    names(CombinedDF) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
		    Simb <- merge(CombinedDF, Xs, by = "Xj")
		 	Simb$HR <- (exp((Simb$Xj - Simb$Xl) * Simb$Coef) - 1) * 100
		    Simb$Comparison <- paste(Simb$Xj, "vs.", Simb$Xl)
		}
	}
	else if (qi == "Hazard Ratio"){
	  	if (length(Xj) != length(Xl)){
	      stop("Xj and Xl must be the same length.")
	    } 
	    else {
		 	Xs <- data.frame(Xj, Xl)   	
			CombinedDF <- MergeX(Xj)
		    names(CombinedDF) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
		    Simb <- merge(CombinedDF, Xs, by = "Xj")
		 	Simb$HR <- exp((Simb$Xj - Simb$Xl) * Simb$Coef)
		    Simb$Comparison <- paste(Simb$Xj, "vs.", Simb$Xl)
	   	}
	}
	else if (qi == "Hazard Rate"){
		Xl <- NULL
		print("Xl is ignored")       
		Simb <- MergeX(Xj)
	    names(Simb) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
	 	Simb$HR <- exp(Simb$Xj * Simb$Coef)
	  	bfit <- basehaz(obj)
	  	## Currently does not support strata
	  	if (!is.null(bfit$strata)){
	  		stop("coxsimSpline currently does not support strata.")
	  	}
	  	bfit$FakeID <- 1
	  	Simb$FakeID <- 1
		bfitDT <- data.table(bfit, key = "FakeID")
		SimbDT <- data.table(Simb, key = "FakeID")
		SimbCombDT <- SimbDT[bfitDT]
		Simb <- data.frame(SimbCombDT)
	  	Simb$HRate <- Simb$hazard * Simb$HR 
	  	Simb <- Simb[, -1]	
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
  class(SimbPerc) <- "simspline"
  SimbPerc
}