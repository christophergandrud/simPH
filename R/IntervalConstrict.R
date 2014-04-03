#' Constrict simulations to a defined interval
#'
#' \code{IntervalConstrict} is an internal function to constrict a set of 
#' simulations to user defined interval.
#'
#' @param Simb character string naming the data frame with the simulations.
#' @param SubVar character vector the variable names to subset the simulations 
#' by.
#' @param qi character vector naming the type of quantity of interest.
#' @param spin logical for whether or not to use the shortest probability 
#' interval or the central interval.
#' @param ci numeric confidence interval measure.
#' @param extremesDrop logical whether or not to drop simulated quantity of interest values that are \code{Inf}, \code{NA}, \code{NaN} and \eqn{> 1000000}. These values are difficult to plot \code{\link{simPH}} and may prevent \code{spin} from finding the central interval.
#' 
#' @importFrom plyr ddply
#' @keywords internals
#' @noRd

IntervalConstrict <- function(Simb = Simb, SubVar = SubVar, qi = qi,
                              spin = FALSE, ci = 0.95, 
                              extremesDrop = extremesDrop)
{
	if (qi == "Hazard Rate" & isTRUE(spin)){
		message(paste("spin currently unsupported for Hazard Rates.",
                      "\nThe central interval will be found instead."))
		spin <- FALSE
	}
    if (isTRUE(extremesDrop)){
        SimbNoExt <- Simb[!is.infinite(Simb$QI), ]
        SimbNoExt <- SimbNoExt[!is.na(SimbNoExt$QI), ]
        SimbNoExt <- SimbNoExt[!is.nan(SimbNoExt$QI), ]
        SimbNoExt <- SimbNoExt[SimbNoExt$QI <= 1000000, ]

        # Messages informing the user of dropped infinite values (if any)
        extDropped <- nrow(Simb) - nrow(SimbNoExt)
        extDroppedPerc <- round((extDropped/nrow(Simb) * 100), digits = 1) 
        if (extDropped > 0){ 
            message(paste0(extDropped, ' (', extDroppedPerc, '%)',
                   ' simulations dropped for having', 
                   '\nquantity of interest values: Inf/NA/NaN/>1000000.\n',
                   '\nTo keep extreme values set extremesDrop = FALSE.\n'))
        }
        Simb <- SimbNoExt
    }
    if (Inf %in% Simb$QI){
        if (isTRUE(spin)){
            stop("spin cannot be TRUE when there are infinite values for your quantity of interest.",
                call. = FALSE)
        } else {
            message(paste("Warning infinite values calculated for your quantity of interest.", 
                "\nConsider changing the difference between Xj and Xl."))
        }
    }
    if (any(Simb$QI > 800) & isTRUE(spin)){
    	message(paste("Warning large quantity of interest values estimated.",
            "\nSPIn may not be found. If so, try spin = FALSE."))
    }

	Lower <- Upper <- NULL
	if (qi == "Relative Hazard" |qi == "Hazard Ratio"| qi == "Hazard Rate"){
		lb <- 0
	} 
	else if (qi == "First Difference"){
		lb <- -100
	} 
	else if (qi == "Marginal Effect"){
		lb <- -Inf
	}

	if (!isTRUE(spin))
	{
		Bottom <- (1 - ci)/2
		Top <- 1 - Bottom
		SimbPerc <- eval(parse(text = 
         paste0("ddply(Simb, SubVar, transform, Lower = QI < quantile(QI,", 
                Bottom, ", na.rm = TRUE))")))
		SimbPerc <- eval(parse(text = 
         paste0("ddply(SimbPerc, SubVar, transform, Upper = QI > quantile(QI,", 
                Top, ", na.rm = TRUE))" )))
	}

	# Drop simulations outside of the shortest probability interval
	else if (isTRUE(spin))
	{
		SimbPerc <- eval(parse(text = paste0("ddply(Simb, SubVar, transform,",
                        "Lower = QI < simPH:::SpinBounds(QI, conf = ", ci,
                         ", lb = ", lb, ", LowUp = 1))" )))
		SimbPerc <- eval(parse(text = paste0("ddply(SimbPerc, SubVar,", 
                        "transform, Upper = QI > simPH:::SpinBounds(QI,",
                        "conf = ", ci, ", lb = ", lb, ", LowUp = 2))" )))
	}

	SimbPerc <- subset(SimbPerc, Lower == FALSE & Upper == FALSE)
	return(SimbPerc)
}

#' Drop any individual simulation with at least one value outside of the
#' specified confidence bounds 
#'
#' @param SimIn data frame
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#'
#' @keywords internals
#' @noRd

OutlierDrop <- function(SimIn){
    # CRAN nonsense
    SimID <- dummy <- Rows <- NULL 
    
    # Drop simulations that do not have all values within
    # the central interval
    SimIn$dummy <- 1
    Sims <- dplyr::group_by(SimIn, SimID) 
    Sims <- dplyr::mutate(Sims, Rows = sum(dummy))
    MaxRows <- max(Sims$Rows)
    Sims <- as.data.frame(Sims)
    Sims <- subset(Sims, Rows == MaxRows)
    Sims
}

