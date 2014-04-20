#' A method for plotting simulation objects created by simPH
#'
#' \code{simGG} a method for ploting simulation objects created by simPH.
#' @param obj an object created by one of simPH's simulation commands.
#' @param ... arguments to be passed to methods.
#'
#' @examples
#' \dontrun{
#'  # Load Carpenter (2002) data
#' data("CarpenterFdaData")
#'
#' # Load survival package
#' library(survival)
#'
#' # Run basic model
#' M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx,
#'             data = CarpenterFdaData)
#'
#' # Simulate Marginal Effect of lethal for multiple
#' # values of prevgenx
#' Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx",
#'                        X2 = seq(2, 115, by = 5), spin = TRUE)
#'
#' # Plot simulations
#' simGG(Sim1)
#' }
#'
#' @seealso \code{\link{simGG.siminteract}}, \code{\link{simGG.simtvc}},
#' \code{\link{simGG.simlinear}}, \code{\link{simGG.simpoly}},
#' \code{\link{simGG.simspline}}
#'
#' @export

simGG <- function(obj, ...){
    UseMethod("simGG", obj)
}
