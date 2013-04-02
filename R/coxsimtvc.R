#' Simulate time-varying hazards from coxph fitted model objects
#' 
#' \code{coxsimtvc} simulates time-varying relative hazards, first differences, and hazard ratios from models estimated with \code{\link{coxph}} using the multivariate normal distribution.
#' @param obj a coxph fitted model object with a time interaction. 
#' @param b the non-time interacted variable's name.
#' @param btvc the time interacted variable's name.
#' @param qi character string indicating what quantity of interest you would like to calculate. Can be \code{'Relative Hazard'}, \code{'First Difference'}, \code{'Hazard Ratio'}, \code{'Hazard Rate'}. Default is \code{qi = 'Relative Hazard'}. If \code{qi = 'First Difference'} or \code{qi = 'Hazard Ratio'} then you can set \code{Xj} and \code{Xl}.
#' @param Xj numeric vector of fitted values for Xj. Must be the same length as Xl. Default is \code{Xj = 1}. Only applies if \code{qi = 'First Difference'} or \code{qi = 'Hazard Ratio'}.
#' @param Xl numeric vector of fitted values for Xl. Must be the same length as Xj. Default is \code{Xl = 0}. Only applies if \code{qi = 'First Difference'} or \code{qi = 'Hazard Ratio'}.
#' @param nsim the number of simulations to run per point in time. Default is \code{nsim = 1000}.
#' @param tfun function of time that btvc was multiplied by. Default is "linear". Can also be "log" (natural log) and "power". If \code{tfun = "power"} then the pow argument needs to be specified also.
#' @param pow if \code{tfun = "power"}, then use pow to specify what power the time interaction was raised to.
#' @param from point in time from when to begin simulating coefficient values
#' @param to point in time to stop simulating coefficient values
#' @param by time intervals by which to simulate coefficient values
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = "95"}, i.e. keep the middle 95 percent. Other possibilities include: \code{"90"}, \code{"99"}, \code{"all"}.
#'
#' @return a simtvc object
#' @details Simulates time-varying relative hazards, first differences, and hazard ratios using parameter estimates from \code{coxph} models. Can also simulate hazard rates for multiple strata.
#'
#' Relative hazards are found using:
#' \deqn{RH = e^{\beta_{1} + \beta_{2}f(t)}}
#' where \eqn{f(t)} is the function of time.
#'
#' First differences are found using:
#' \deqn{FD = (e^{(X_{j} - X_{l}) (\beta_{1} + \beta_{2}f(t))} - 1) * 100}
#' where \eqn{X_{j}} and \eqn{X_{l}} are some values of \eqn{X} to contrast.
#'
#' Hazard ratios are calculated using:
#' \deqn{FD = e^{(X_{j} - X_{l}) (\beta_{1} + \beta_{2}f(t))}}
#' When simulating non-stratifed time-varying harzards \code{coxsimtvc} uses the point estimates for a given coefficient \eqn{\hat{\beta}_{x}} and its time interaction \eqn{\hat{\beta}_{xt}} along with the variance matrix (\eqn{\hat{V}(\hat{\beta})}) estimated from a \code{coxph} model. These are used to draw values of \eqn{\beta_{x}} and \eqn{\beta_{xt}} from the multivariate normal distribution \eqn{N(\hat{\beta},\: \hat{V}(\hat{beta}))}.
#'
#' When simulating stratified time-varying hazard rates \eqn{H} for a given strata \eqn{k}, \code{coxsimtvc} uses:
#' \deqn{H_{kxt} = \hat{\beta_{k0t}}\exp{\hat{\beta_{x}} + \beta_{xt}(t)}}
#' The resulting simulation values can be plotted using \code{\link{ggtvc}}. 
#'
#' @examples
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#' 
#' # Load survival package
#' library(survival)
#' 
#' # Create natural log time interactions
#' Golubtvc <- function(x){
#'   assign(paste0("l", x), tvc(GolubEUPData, b = x, tvar = "end", tfun = "log"))
#' }
#' 
#' GolubEUPData$Lcoop <-Golubtvc("coop")
#' GolubEUPData$Lqmv <- Golubtvc("qmv")
#' GolubEUPData$Lbacklog <- Golubtvc("backlog")
#' GolubEUPData$Lcodec <- Golubtvc("codec")
#' GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
#' GolubEUPData$Lthatcher <- Golubtvc("thatcher") 
#' 
#' # Run Cox PH Model
#' M1 <- coxph(Surv(begin, end, event) ~ 
#'             qmv + qmvpostsea + qmvpostteu + 
#'             coop + codec + eu9 + eu10 + eu12 +
#'             eu15 + thatcher + agenda + backlog +
#'             Lqmv + Lqmvpostsea + Lcoop + Lcodec +
#'             Lthatcher + Lbacklog, 
#'          data = GolubEUPData,
#'          ties = "efron")
#'          
#' # Create simtvc object for Relative Hazard
#' Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
#'                    tfun = "log", from = 80, to = 2000, 
#'                    by = 15, ci = "99")
#' 
#' # Create simtvc object for First Difference  
#' Sim2 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
#'                   qi = "First Difference", 
#'                   tfun = "log", from = 80, to = 2000, 
#'                   by = 15, ci = "99")
#' 
#' # Create simtvc object for Hazard Ratio  
#' Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
#'                   qi = "Hazard Ratio", Xj = c(191, 229), 
#'                   Xl = c(0, 0),
#'                   tfun = "log", from = 80, to = 2000, 
#'                   by = 15, ci = "99")
#'
#' @seealso \code{\link{ggtvc}}, \code{\link{rmultinorm}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#' @import MSBVAR plyr reshape2 survival data.table
#' @export
#' @references Golub, Jonathan, and Bernard Steunenberg. 2007. “How Time Affects EU Decision-Making.” European Union Politics 8(4): 555–66.
#'
#' Licht, Amanda A. 2011. “Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.” Political Analysis 19: 227–43.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. “Making the Most of Statistical Analyses: Improving Interpretation and Presentation.” American Journal of Political Science 44(2): 347–61.

coxsimtvc <- function(obj, b, btvc, qi = "Relative Hazard", Xj = 1, Xl = 0, tfun = "linear", pow = NULL, nsim = 1000, from, to, by, ci = "95")
{
  # Parameter estimates & Varance/Covariance matrix
  Coef <- matrix(obj$coefficients)
  VC <- vcov(obj)
    
  # Draw values from the multivariate normal distribution
  Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
  DrawnDF <- data.frame(Drawn)
 
  # Extract simulations for variables of interest
  dfn <- names(DrawnDF)
  bpos <- match(b, dfn)
  btvcpos <- match(btvc, dfn)
  
  Drawn <- data.frame(Drawn[, c(bpos, btvcpos)])
  Drawn$ID <- 1:nsim
  
  # Create time function
  tfunOpts <- c("linear", "log", "power")
  TestforTOpts <- tfun %in% tfunOpts
  if (TestforTOpts == FALSE){
    stop("Must specify tfun as 'linear', 'log', or 'power'")
  }
    
  if (tfun == "linear"){
    tf <- seq(from = from, to = to, by = by)
  } else if (tfun == "log"){
    tf <- log(seq(from = from, to = to, by = by))
  } else if (tfun == "power"){
    tf <- (seq(from = from, to = to, by = by))^pow
  }
  
  # Multiply time function with btvc
  TVSim <- outer(Drawn[,2], tf)
  TVSim <- data.frame(melt(TVSim))
  names(TVSim) <- c("ID", "time", "TVC")
  time <- 1:length(tf)
  TempDF <- data.frame(time, tf)
  TVSim <- merge(TVSim, TempDF)

  # Combine with non TVC version of the variable
  TVSim <- merge(Drawn, TVSim, by = "ID")
  TVSim$CombCoef <- TVSim[[2]] + TVSim$TVC

  # Find quantity of interest
  if (qi == "Relative Hazard"){
      message("All Xl ignored.")
      Xs <- data.frame(Xj)
      names(Xs) <- c("Xj")
      Xs$Comparison <- paste(Xs[, 1])
      TVSim <- merge(TVSim, Xs)
      TVSim$HR <- exp(TVSim$CombCoef * TVSim$Xj)
  } else if (qi == "First Difference"){
    if (length(Xj) != length(Xl)){
      stop("Xj and Xl must be the same length.")
    } else {
      TVSim$HR <- exp(TVSim$CombCoef)
      Xs <- data.frame(Xj, Xl)
      Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
      TVSim <- merge(TVSim, Xs)
      TVSim$FirstDiff <- (exp((TVSim$Xj - TVSim$Xl) * TVSim$CombCoef) - 1) * 100
    }
  } else if (qi == "Hazard Ratio"){
    if (length(Xl) > 1 & length(Xj) != length(Xl)){
      stop("Xj and Xl must be the same length.")
    }
    else {
      if (length(Xj) > 1 & length(Xl) == 1){
        Xl <- rep(0, length(Xj))
      }
      Xs <- data.frame(Xj, Xl)
      Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
      TVSim <- merge(TVSim, Xs)
      TVSim$HR <- exp((TVSim$Xj - TVSim$Xl) * TVSim$CombCoef)
    }
  } else if (qi == "Hazard Rate"){
      Xl <- NULL
      message("Xl is ignored") 
      Xs <- data.frame(Xj)
      Xs$HRValue <- paste(Xs[, 1])
      Simb <- merge(TVSim, Xs)
      Simb$HR <- exp(Simb$Xj * Simb$CombCoef)  
      bfit <- basehaz(obj)
      bfit$FakeID <- 1
      Simb$FakeID <- 1
      bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
      SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
      SimbCombDT <- SimbDT[bfitDT, allow.cartesian = TRUE]
      Simb <- data.frame(SimbCombDT)
      Simb$HRate <- Simb$hazard * Simb$HR 
      TVSim <- Simb[, -1]
  }

  # Drop simulations outside of 'confidence bounds'
  TVSim <- TVSim[order(TVSim$time),]
  if (ci == "all"){
    TVSimPerc <- TVSim 
  } else if (ci == "90"){
    TVSimPerc <- ddply(TVSim, .(time, Xj), transform, Lower = HR < quantile(HR, c(0.05)))
    TVSimPerc <- ddply(TVSimPerc, .(time, Xj), transform, Upper = HR > quantile(HR, c(0.95)))
  } else if (ci == "95"){
    TVSimPerc <- ddply(TVSim, .(time, Xj), transform, Lower = HR < quantile(HR, c(0.025)))
    TVSimPerc <- ddply(TVSimPerc, .(time, Xj), transform, Upper = HR > quantile(HR, c(0.975)))
  } else if (ci == "99"){
    TVSimPerc <- ddply(TVSim, .(time, Xj), transform, Lower = HR < quantile(HR, c(0.005)))
    TVSimPerc <- ddply(TVSimPerc, .(time, Xj), transform, Upper = HR > quantile(HR, c(0.995)))
  }
  if (ci != "all"){
    TVSimPerc <- subset(TVSimPerc, Lower == FALSE & Upper == FALSE)
  }

  # Create real time variable
  if (tfun == "linear"){
    TVSimPerc$RealTime <- TVSimPerc$tf
  } else if (tfun == "log"){
    TVSimPerc$RealTime <- exp(TVSimPerc$tf)
  } else if (tfun == "power"){
    TVSimPerc$RealTime <- TVSimPerc$tf^(1/pow)
  }
  
  # Final clean up
  class(TVSimPerc) <- "simtvc"
  TVSimPerc
}