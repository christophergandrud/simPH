#' A method for plotting simulation objects created by simPH.
#'
#' @param obj an object created by one of simPH's simulation commands.
#' @param ... arguments to be passed to methods.
#' 
#' @seealso \code{\link{simGG.siminteract}}, \code{\link{simGG.simtvc}}, \code{\link{simGG.simlinear}}, \code{\link{simGG.simpoly}}, \code{\link{simGG.simspline}}
#' @import ggplot2
#' @export simGG  
#' 
simGG <- function(obj, ...){
  UseMethod("simGG", obj)
}