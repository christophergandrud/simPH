#' Transform the simulation object to include only the min and max of the
#' constricted interval, as well as the lower and upper bounds of the middle 50
#' percent of the constricted interval
#'
#' \code{MinMaxLines} is an internal function to transform the simulation
#' object to include only the min and max of the interval set by \code{ci} in
#' the \code{coxsim} command, as well as the lower and upper bounds of the
#' middle 50 percent of this interval. It also returns the median.
#'
#' @param df a data frame or a simulation class object.
#' @param byVars character vector of the variables to subset the data frame by.
#' The default is \code{'Xj'}.
#' @param hr logical indicating whether or not \code{df} contains a hazard rate.
#' @param strata logical indicating whether or not \code{df} contains a
#' stratified hazard rate.
#' @param clean logical, whether or not to clean up the output data frame to
#' only include \code{byVars}, \code{Min_CI}, \code{Lower50_CI}, \code{median},
#' \code{Upper50_CI}, \code{Max_CI}.
#'
#' @importFrom plyr ddply
#' @keywords internals
#' @export

MinMaxLines <- function(df, byVars = "Xj", hr = FALSE, strata = FALSE,
                        clean = FALSE){
    Xj <- QI <- Time <- HRValue <- HRate <- Strata <- NULL
	if (class(df) != 'data.frame') class(df) <- 'data.frame'
    if (isTRUE(hr) & !isTRUE(strata)){
        byVars <- c("Time", "HRValue")
    }
    else if (isTRUE(hr) & !isTRUE(strata)){
        byVars <- c("Time", "HRValue", "Strata")
    }
    if (!isTRUE(hr)){
        Linesdf <- ddply(df, byVars, transform, Median = median(QI))
        Linesdf <- ddply(Linesdf, byVars, transform, Max = max(QI))
        Linesdf <- ddply(Linesdf, byVars, transform, Min = min(QI))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Lower50 = quantile(QI, 0.25))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Upper50 = quantile(QI, 0.75))

        Linesdf <- Linesdf[!duplicated(Linesdf[, byVars]), ]
    }
    else if (isTRUE(hr) & !isTRUE(strata)){
        Linesdf <- ddply(df, byVars, transform, Median = median(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform, Max = max(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform, Min = min(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Lower50 = quantile(HRate, 0.25))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Upper50 = quantile(HRate, 0.75))

        Linesdf <- Linesdf[!duplicated(Linesdf[, c(1, 3)]), ]

    }
    else if (isTRUE(hr) & isTRUE(strata)){
        Linesdf <- ddply(df, byVars, transform, Median = median(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform, Max = max(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform, Min = min(HRate))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Lower50 = quantile(HRate, 0.25))
        Linesdf <- ddply(Linesdf, byVars, transform,
                        Upper50 = quantile(HRate, 0.75))

        Linesdf <- Linesdf[!duplicated(
                            Linesdf[, c("Time", "HRValue", "Strata")]), ]
    }
    if (isTRUE(clean)){
        Linesdf <- Linesdf[, c(byVars, 'Min', 'Lower50', 'Median', 'Upper50',
                                'Max')]
        names(Linesdf) <- c(byVars, 'Min_CI', 'Lower50_CI','Median',
                            'Upper50_CI', 'Max_CI')
        }
    return(Linesdf)
}


#' Create a variable of each simulation value's percentile in the distribution.
#'
#' @param SimIn data frame of simulations.
#' @param xaxis character string. The column that will form the x-axis in
#' the plot.
#' @param yaxis character string. The column that will form the y-axis in
#' the plot.
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#'
#' @keywords internals
#' @noRd

PercRank <- function(SimIn, xaxis, yaxis = 'QI'){
    Xaxis <- QI <- NULL

    names(SimIn)[names(SimIn) == xaxis] <- 'Xaxis'
    names(SimIn)[names(SimIn) == yaxis] <- 'QI'

    PlainPercRank <- function(x) trunc(rank(x))/length(x)

    Temp <- dplyr::group_by(SimIn, Xaxis)
    Temp <- dplyr::mutate(Temp, PercRank = PlainPercRank(QI))

    # Center on the 50th percentile
    Temp$PercRank <- abs(Temp$PercRank - 0.5)
    Temp$PercRank <- round(abs(0.5 - Temp$PercRank), 1)

    names(Temp)[names(Temp) == 'Xaxis'] <- xaxis
    names(Temp)[names(Temp) == 'QI'] <- yaxis
    Temp
}
