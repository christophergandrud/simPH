#' Create a time interaction variable
#' 
#' \code{tvc} creates a time interaction variable that can be used in a coxph model (or any other model with time interactions)
#' @param data a data frame
#' @param b the non-time interacted variable's name
#' @param tvar the time variable's name
#' @param tfun function of time that btvc was multiplied by. Default is \code{tfun = "linear"}. Can also be \code{tfun = "log"} (natural log) and \code{tfun = "power"}. If \code{tfun = "power"} then the pow argument needs to be specified also.
#' @param pow if \code{tfun = "power"}, then use pow to specify what power to raise the time interaction to.
#' @return a vector
#' @details Interacts a variable with a specified function of time. Possible functions of time include \code{'linear'}, natural \code{'log'}, and exponentiated (\code{'power'}).
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#' 
#' # Create natural log time interaction with the qmv variable
#' GolubEUPData$Lqmv <- tvc(GolubEUPData, b = "qmv", tvar = "end", tfun = "log")
#' @seealso \code{\link{ggtvc}}, \code{coxsimtvc}, \code{\link{survival}}, and \code{\link{coxph}}
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
