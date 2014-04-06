#' Create a time interaction variable
#' 
#' \code{tvc} creates a time interaction variable that can be used in a coxph 
#' model (or any other model with time interactions)
#' @param data a data frame
#' @param b the non-time interacted variable's name. Either a single value or a 
#' vector of variable names can be entered.
#' @param tvar the time variable's name
#' @param tfun function of time that btvc was multiplied by. Default is 
#' \code{tfun = "linear"}. Can also be \code{tfun = 'log'} (natural log) and 
#' \code{tfun = 'power'}. If \code{tfun = 'power'} then the pow argument needs 
#' to be specified also.
#' @param pow if \code{tfun = 'power'}, then use pow to specify what power to 
#' raise the time interaction to.
#' @param vector logical. Whether or not to return one vector a or a data frame.
#' Can only be used if only one \code{b} is included.
#' @return a data frame or vector. If a data frame is returned it will include 
#' all of the original variables as well as the interactions denoted by a 
#' variable name '\code{bn}_\code{tfun}', where \code{bn} is one variable name 
#' from \code{b} and \code{tfun} as entered into the function.
#' @details Interacts a variable with a specified function of time. Possible 
#' functions of time include \code{'linear'}, natural \code{'log'}, and 
#' exponentiated (\code{'power'}).
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data('GolubEUPData')
#' 
#' # Subset PURELY TO SPEED UP THE EXAMPLE
#' GolubEUPData <- GolubEUPData[1:500, ]
#'
#' # Expand data into equally spaced time intervals
#' GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
#'                   Time = 'begin', Time2 = 'end', event = 'event') 
#' 
#' # Create natural log time interaction with the qmv variable
#' GolubEUPData$Lqmv <- tvc(GolubEUPData, b = 'qmv', tvar = 'end', tfun = 'log', 
#'                          vector = TRUE)
#' 
#' # Create interactions for a vector of variables
#' BaseVars <- c('qmv', 'backlog', 'coop', 'codec', 'qmvpostsea', 'thatcher')
#' Test <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')
#'                                           
#' @seealso \code{\link{SurvExpand}}, \code{\link{simGG.simtvc}}, 
#' \code{\link{coxsimtvc}}, \code{\link{survival}}, and \code{\link{coxph}}
#' @keywords utilities
#' @export


tvc <- function(data, b, tvar, tfun = "linear", pow = NULL, vector = FALSE)
{
    NumBs <- length(b)
    # Errors
    if (NumBs != 1 & isTRUE(vector)){
      warning('If b > 1 then vector is reset to FALSE.', call. = FALSE)
      vector <- FALSE
    }

    tfunOpts <- c("linear", "log", "power")
    TestforTOpts <- tfun %in% tfunOpts
    if (!isTRUE(TestforTOpts)){
      stop("Must specify tfun as 'linear', 'log', or 'power'.", , call. = FALSE)
    }

    # Function to create interactions
    tvcCreate <- function(data, b, tvar, tfun){
        if (tfun == "linear"){
            data[, b] * data[, tvar]
        } else if (tfun == "log"){
            data[, b] * log(data[, tvar])
        } else if (tfun == "power") {
            data[, b] * (data[, tvar])^pow
        }
    }  

    if (isTRUE(vector)){
       Out <- tvcCreate(data = data, b = b, tvar = tvar, tfun = tfun)
    }
    else if (!isTRUE(vector)){
      for (i in b){
        newVar <- paste0(i, "_", tfun)
        data[, newVar] <- tvcCreate(data = data, b = i, tvar = tvar, tfun = tfun)
      }
      Out <- data
    }
    return(Out)
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