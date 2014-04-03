#' Create a time interaction variable
#' 
#' \code{tvc} creates a time interaction variable that can be used in a coxph 
#' model (or any other model with time interactions)
#' @param data a data frame
#' @param b the non-time interacted variable's name
#' @param tvar the time variable's name
#' @param tfun function of time that btvc was multiplied by. Default is 
#' \code{tfun = "linear"}. Can also be \code{tfun = "log"} (natural log) and 
#' \code{tfun = "power"}. If \code{tfun = "power"} then the pow argument needs 
#' to be specified also.
#' @param pow if \code{tfun = "power"}, then use pow to specify what power to 
#' raise the time interaction to.
#' @return a vector
#' @details Interacts a variable with a specified function of time. Possible 
#' functions of time include \code{'linear'}, natural \code{'log'}, and 
#' exponentiated (\code{'power'}).
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#' 
#' # Subset PURELY TO SPEED UP THE EXAMPLE
#' GolubEUPData <- GolubEUPData[1:500, ]
#'
#' # Expand data into equally spaced time intervals
#' GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
#'                   Time = 'begin', Time2 = 'end', event = 'event') 
#' 
#' # Create natural log time interaction with the qmv variable
#' GolubEUPData$Lqmv <- tvc(GolubEUPData, b = "qmv", 
#'                        tvar = "end", tfun = "log")
#' @seealso \code{\link{SurvExpand}}, \code{\link{simGG.simtvc}}, 
#' \code{coxsimtvc}, \code{\link{survival}}, and \code{\link{coxph}}
#' @keywords utilities
#' @export


tvc <- function(data, b, tvar, tfun = "linear", pow = NULL)
{
  dfn <- names(data)
  bpos <- match(b, dfn)
  tvarpos <- match(tvar, dfn)
  
  tfunOpts <- c("linear", "log", "power")
  TestforTOpts <- tfun %in% tfunOpts
  if (!isTRUE(TestforTOpts)){
    stop("Must specify tfun as 'linear', 'log', or 'power'")
  }
    
  if (tfun == "linear"){
    data[[bpos]] * data[[tvarpos]]
  } else if (tfun == "log"){
    data[[b]] * log(data[[tvarpos]])
  } else if (tfun == "power") {
    data[[b]] * (data[[tvarpos]])^pow
  }
}

#' Create a sequence of Xl values
#' 
#' \code{setXl} creates a sequence of \code{Xl} values given a sequence of 
#' \code{Xj} values and a fixed difference.
#' @param Xj numeric vector of fitted values for the covariate of interest to 
#' simulate for.
#' @param diff numeric vector of length 1. It specifies the difference between 
#' \code{Xj} and \code{Xl}. \code{Xl} is always smaller than \code{Xj}.
#' 
#' @return a vector
#' 
#' @examples 
#' # Set Xj
#' setXj = seq(1100, 1700, by = 10)
#' 
#' # Find Xl that are 1 less than Xj 
#' setXl(Xj = setXj, diff = 1) 
#' @keywords utilities
#' @export

setXl <- function(Xj, diff = 1){
  # Errors
  if (class(Xj) != 'numeric') stop('Xj must be numeric',
                                                call. = FALSE)
  if (class(diff) != 'numeric') stop('diff must be numeric',
                                     call. = FALSE)
  if (length(diff) != 1) stop('diff can only be one value', call. = FALSE)
  if (diff <= 0){
    message('diff cannot be negative. Using the absolute value of diff.', 
      call. = FALSE)
    diff <- abs(diff)
  } 
  if (diff == 0) stop('diff cannot be 0.', call. = FALSE)

  # Set Xl  
  set <- Xj - diff
  return(set)
}