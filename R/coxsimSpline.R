#' Simulate quantities of interest for penalized splines from Cox Proportional Hazards models
#'
#' \code{coxsimSpline} simulates quantities of interest from penalized splines using multivariate normal distributions.
#' @param obj a \code{\link{coxph}} class fitted model object with a penalized spline. These can be plotted with \code{\link{simGG}}.
#' @param bspline a character string of the full \code{\link{pspline}} call used in \code{obj}. It should be exactly the same as how you entered it in \code{\link{coxph}}. You also need to enter a white spece before and after all equal (\code{=}) signs.
#' @param bdata a numeric vector of splined variable's values.
#' @param qi quantity of interest to simulate. Values can be \code{"Relative Hazard"}, \code{"First Difference"}, \code{"Hazard Ratio"}, and \code{"Hazard Rate"}. The default is \code{qi = "Relative Hazard"}. Think carefully before using \code{qi = "Hazard Rate"}. You may be creating very many simulated values which can be very computationally intensive to do. Adjust the number of simulations per fitted value with \code{nsim}.
#' @param Xj numeric vector of fitted values for \code{b} to simulate for.
#' @param Xl numeric vector of values to compare \code{Xj} to. Note if \code{qi = "Relative Hazard"} or \code{"Hazard Rate"} only \code{Xj} is relevant.
#' @param nsim the number of simulations to run per value of \code{Xj}. Default is \code{nsim = 1000}.
#' @param ci the proportion of simulations to keep. The default is \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE} then \code{ci} is the confidence level of the shortest probability interval. Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest probability interval rather than the middle simulations.
#'
#' @return a \code{simspline} object
#'
#' @details Simulates relative hazards, first differences, hazard ratios, and hazard rates for penalized splines from Cox Proportional Hazards models. These can be plotted with \code{\link{simGG}}. 
#' A Cox PH model with one penalized spline is given by:
#'  \deqn{h(t|\mathbf{X}_{i})=h_{0}(t)\mathrm{e}^{g(x)}}
#'
#' where \eqn{g(x)} is the penalized spline function. For our post-estimation purposes \eqn{g(x)} is basically a series of linearly combined coefficients such that:
#'
#'  \deqn{g(x) = \beta_{k_{1}}(x)_{1+} + \beta_{k_{2}}(x)_{2+} + \beta_{k_{3}}(x)_{3+} + \ldots + \beta_{k_{n}}(x)_{n+}}
#'
#' where \eqn{k} are the equally spaced spline knots with values inside of the range of observed \eqn{x} and \eqn{n} is the number of knots. 
#'
#' We can again draw values of each \eqn{\beta_{k_{1}}, \ldots \beta_{k_{n}}} from the multivariate normal distribution described above. We then use these simulated coefficients to estimates quantities of interest for a range covariate values. For example, the first difference between two values \eqn{x_{j}} and \eqn{x_{l}} is:
#'
#'    \deqn{\%\triangle h_{i}(t) = (\mathrm{e}^{g(x_{j}) - g(x_{l})} - 1) * 100}
#'
#' Relative hazards and hazard ratios can be calculated by extension.
#'
#' Currently \code{coxsimSpline} does not support simulating hazard rates form multiple stratified models.
#'
#' @examples
#' ## dontrun
#' # Load Carpenter (2002) data
#' # data("CarpenterFdaData")
#' 
#' # Load survival package
#' # library(survival)
#' 
#' # Run basic model
#' # From Keele (2010) replication data
#' # M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 + 
#' #				acutediz + hosp01  + pspline(hospdisc, df = 4) + 
#' #				pspline(hhosleng, df = 4) + mandiz01 + femdiz01 + peddiz01 +
#' #				orphdum + natreg + vandavg3 + wpnoavg3 + 
#' #				pspline(condavg3, df = 4) + pspline(orderent, df = 4) + 
#'#				pspline(stafcder, df = 4), data = CarpenterFdaData)
#'
#' # Simulate Relative Hazards for orderent
#' # Sim1 <- coxsimSpline(M1, bspline = "pspline(stafcder, df = 4)", 
#' #                    bdata = CarpenterFdaData$stafcder,
#' #                    qi = "Hazard Ratio",
#' #                    Xj = seq(1100, 1700, by = 10), 
#' #                    Xl = seq(1099, 1699, by = 10), spin = TRUE)
#'
#' ## dontrun  
#' # Simulate Hazard Rates for orderent
#' # Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#' #                       bdata = CarpenterFdaData$orderent,
#' #                       qi = "Hazard Rate",
#' #                       Xj = seq(2, 53, by = 3),
#' #                       nsim = 100)
#'
#'
#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#'
#' @references Luke Keele, "Replication data for: Proportionally Difficult: Testing for Nonproportional Hazards In Cox Models", 2010, \url{http://hdl.handle.net/1902.1/17068} V1 [Version]. 
#' 
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' American Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient Shortest Probability Intervals.'' Arvix. \url{http://arxiv.org/pdf/1302.2142v1.pdf}.
#' 
#' @import data.table
#' @importFrom stringr word str_match str_replace
#' @importFrom reshape2 melt
#' @importFrom survival basehaz
#' @importFrom MSBVAR rmultnorm
#' @export

coxsimSpline <- function(obj, bspline, bdata, qi = "Relative Hazard", Xj = 1, Xl = 0, nsim = 1000, ci = 0.95, spin = FALSE)
{ 
	QI <- NULL
	# Ensure that qi is valid
	qiOpts <- c("Relative Hazard", "First Difference", "Hazard Rate", "Hazard Ratio")
	TestqiOpts <- qi %in% qiOpts
	if (!isTRUE(TestqiOpts)){
		stop("Invalid qi type. qi must be 'Relative Hazard', 'Hazard Rate', 'First Difference', or 'Hazard Ratio'")
	}

	if (nsim > 10 & qi == "Hazard Rate"){
		message(paste0("Warning: finding Hazard Rates with ", nsim, " simulations may take awhile.  Consider decreasing the number of simulations with nsim."))
	}

	if (is.null(Xl) & qi != "Hazard Rate"){
		Xl <- rep(0, length(Xj))
		message("All Xl set to 0.")
	} else if (!is.null(Xl) & qi == "Relative Hazard") {
		Xl <- rep(0, length(Xj))
		message("All Xl set to 0.")
	}

	# Standardise pspline term with white space around '='
	NoSpace1 <- grepl("[^ ]=", bspline)
	NoSpace2 <- grepl("[^ ]=", bspline)
	if (any(NoSpace1) | any(NoSpace2)){
		stop(paste0("Place a white space before and after equal (=) signs in ", bspline, "."))
	}


	# Find term number
	TermNum <- names(obj$pterms)
	bterm <- match(bspline, TermNum)
	if (is.na(bterm)){
		stop(paste0("Unable to find ", bspline, "."))
	}

	# Extract boundary knots for default Boundary.knots = range(x) & number of knots
	#### Note: these can also be found with get("cbase", environment(obj$printfun[[1]])) # (replace 1 with the spline term number
	#### From: http://r.789695.n4.nabble.com/help-on-pspline-in-coxph-td3431829.html 
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
	bword <- word(bspline, 1)
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
		X <- NULL
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
	if (qi == "Relative Hazard" | qi == "Hazard Ratio"){
	  	if (length(Xj) != length(Xl)){
	      stop("Xj and Xl must be the same length.")
	    } 

		Simbj <- MergeX(Xj)
	    names(Simbj) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
		Simbl <- MergeX(Xl)
	    names(Simbl) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xl")

		if (qi == "Hazard Ratio"){
		 	Xs <- data.frame(Xj, Xl)   	
		    Simbj <- merge(Simbj, Xs, by = "Xj")
		    Simbj$Comparison <- paste(Simbj$Xj, "vs.", Simbj$Xl)
		}

	    Simbj$QI <- exp((Simbj$Xj * Simbj$Coef) - (Simbl$Xl * Simbl$Coef))
	    Simb <- Simbj
	}
	else if (qi == "First Difference"){
	  	if (length(Xj) != length(Xl)){
	      stop("Xj and Xl must be the same length.")
	    } 
	    else {
			Simbj <- MergeX(Xj)
		    names(Simbj) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xj")
			Simbl <- MergeX(Xl)
		    names(Simbl) <- c("CoefName", "Coef", "IntervalStart", "IntervalFinish", "Xl")
		 	
		 	Xs <- data.frame(Xj, Xl)   	
		    Simbj <- merge(Simbj, Xs, by = "Xj")
		    Simbj$Comparison <- paste(Simbj$Xj, "vs.", Simbj$Xl)

		 	Simbj$QI <- (exp((Simbj$Xj * Simbj$Coef) - (Simbl$Xl * Simbl$Coef)) - 1) * 100
			Simb <- Simbj
		}
	}
	else if (qi == "Hazard Rate"){
		Xl <- NULL
      message("Xl is ignored. All variables' values other than b fitted at 0.") 
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
		bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
		SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
		SimbCombDT <- SimbDT[bfitDT, allow.cartesian = TRUE]
		Simb <- data.frame(SimbCombDT)
	  	Simb$QI <- Simb$hazard * Simb$HR 
	  	Simb <- Simb[, -1]	
	}

	# Drop simulations outside of 'confidence bounds'
	if (qi != "Hazard Rate"){
		SubVar <- "Xj"
	} else if (qi == "Hazard Rate"){
		SubVar <- c("time", "Xj")
	}

  if (Inf %in% Simb$QI){
  	if (isTRUE(spin)){
		stop("spin cannot be TRUE when there are infinite values for your quantitiy of interest.")
  	} else {
	  	message("Warning infinite values calculated for your quantity of interest. Consider changing the difference between Xj and Xl.")
  	}
  }

  SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar, qi = qi,
          QI = QI, spin = spin, ci = ci)	


  # Final clean up
  class(SimbPerc) <- c("simspline", qi)
  SimbPerc
}
