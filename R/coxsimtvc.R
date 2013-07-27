#' Simulate time-varying quantities of interest from coxph fitted model objects
#' 
#' \code{coxsimtvc} simulates time-varying relative hazards, first differences, and hazard ratios from models estimated with \code{\link{coxph}} using the multivariate normal distribution.
#' @param obj a coxph fitted model object with a time interaction. 
#' @param b the non-time interacted variable's name.
#' @param btvc the time interacted variable's name.
#' @param qi character string indicating what quantity of interest you would like to calculate. Can be \code{'Relative Hazard'}, \code{'First Difference'}, \code{'Hazard Ratio'}, \code{'Hazard Rate'}. Default is \code{qi = 'Relative Hazard'}. If \code{qi = 'First Difference'} or \code{qi = 'Hazard Ratio'} then you can set \code{Xj} and \code{Xl}.
#' @param Xj numeric vector of fitted values for Xj. Must be the same length as \code{Xl} or \code{Xl} must be \code{NULL}. 
#' @param Xl numeric vector of fitted values for Xl. Must be the same length as Xj. Only applies if \code{qi = 'First Difference'} or \code{qi = 'Hazard Ratio'}.
#' @param nsim the number of simulations to run per point in time. Default is \code{nsim = 1000}.
#' @param tfun function of time that btvc was multiplied by. Default is "linear". Can also be "log" (natural log) and "power". If \code{tfun = "power"} then the pow argument needs to be specified also.
#' @param pow if \code{tfun = "power"}, then use pow to specify what power the time interaction was raised to.
#' @param from point in time from when to begin simulating coefficient values
#' @param to point in time to stop simulating coefficient values
#' @param by time intervals by which to simulate coefficient values
#' @param ci the proportion of middle simulations to keep. The default is \code{ci = 0.95}, i.e. keep the middle 95 percent. If \code{spin = TRUE} then \code{ci} is the convidence level of the shortest probability interval. Any value from 0 through 1 may be used.
#' @param spin logical, whether or not to keep only the shortest proability interval rather than the middle simulations.
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
#' When simulating non-stratifed time-varying harzards \code{coxsimtvc} uses the point estimates for a given coefficient \eqn{\hat{\beta}_{x}} and its time interaction \eqn{\hat{\beta}_{xt}} along with the variance matrix (\eqn{\hat{V}(\hat{\beta})}) estimated from a \code{coxph} model. These are used to draw values of \eqn{\beta_{1}} and \eqn{\beta_{2}} from the multivariate normal distribution \eqn{N(\hat{\beta},\: \hat{V}(\hat{beta}))}.
#'
#' When simulating stratified time-varying hazard rates \eqn{H} for a given strata \eqn{k}, \code{coxsimtvc} uses:
#' \deqn{H_{kxt} = \hat{\beta_{k0t}}\exp{\hat{\beta_{1}} + \beta_{2}f(t)}}
#' The resulting simulation values can be plotted using \code{\link{simGG}}. 
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
#'                    Xj = 1, by = 15, ci = 0.99)
#' 
#' # Create simtvc object for First Difference  
#'Sim2 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
#'                  qi = "First Difference", Xj = 1,
#'                  tfun = "log", from = 80, to = 2000,
#'                  by = 15, ci = 0.95)
#' 
#' # Create simtvc object for Hazard Ratio  
#' Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
#'                   qi = "Hazard Ratio", Xj = c(191, 229), 
#'                   Xl = c(0, 0),
#'                   tfun = "log", from = 80, to = 2000, 
#'                   by = 15, ci = 0.5)
#'

#' @seealso \code{\link{simGG}}, \code{\link{survival}}, \code{\link{strata}}, and \code{\link{coxph}}
#'
#' @import data.table
#' @importFrom reshape2 melt
#' @importFrom survival basehaz
#' @importFrom MSBVAR rmultnorm
#' @export
#'
#' @references Golub, Jonathan, and Bernard Steunenberg. 2007. ''How Time Affects EU Decision-Making.'' European Union Politics 8(4): 555-66.
#'
#' Licht, Amanda A. 2011. ''Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.'' Political Analysis 19: 227-43.
#'
#' King, Gary, Michael Tomz, and Jason Wittenberg. 2000. ''Making the Most of Statistical Analyses: Improving Interpretation and Presentation.'' American Journal of Political Science 44(2): 347-61.
#'
#' Liu, Ying, Andrew Gelman, and Tian Zheng. 2013. ''Simulation-Efficient Shortest Probablility Intervals.'' Arvix. http://arxiv.org/pdf/1302.2142v1.pdf.


coxsimtvc <- function(obj, b, btvc, qi = "Relative Hazard", Xj = NULL, Xl = NULL, tfun = "linear", pow = NULL, nsim = 1000, from, to, by, ci = 0.95, spin = FALSE)
{
  QI <- NULL
  ############ Means Not Supported Yet for coxsimtvc ########
  #### means code is a place holder for future versions #####
  means <- FALSE
  # Ensure that qi is valid
  qiOpts <- c("Relative Hazard", "First Difference", "Hazard Rate", "Hazard Ratio")
  TestqiOpts <- qi %in% qiOpts
  if (!isTRUE(TestqiOpts)){
    stop("Invalid qi type. qi must be 'Relative Hazard', 'First Difference', 'Hazard Rate', or 'Hazard Ratio'")
  }

  MeansMessage <- NULL
  if (isTRUE(means) & length(obj$coefficients) == 3){
    means <- FALSE
    MeansMessage <- FALSE
    message("Note: means reset to FALSE. The model only includes the interaction variables.")
  } else if (isTRUE(means) & length(obj$coefficients) > 3){
    MeansMessage <- TRUE
  }

  if (is.null(Xl) & qi != "Hazard Rate"){
    Xl <- rep(0, length(Xj))
    message("All Xl set to 0.")
  } else if (!is.null(Xl) & qi == "Relative Hazard") {
    message("All Xl set to 0.")
  }

  # Create time function
  tfunOpts <- c("linear", "log", "power")
  TestforTOpts <- tfun %in% tfunOpts
  if (!isTRUE(TestforTOpts)){
    stop("Must specify tfun as 'linear', 'log', or 'power'")
  }
    
  if (tfun == "linear"){
    tf <- seq(from = from, to = to, by = by)
  } else if (tfun == "log"){
    tf <- log(seq(from = from, to = to, by = by))
  } else if (tfun == "power"){
    tf <- (seq(from = from, to = to, by = by))^pow
  }

  # Parameter estimates & Varance/Covariance matrix
  Coef <- matrix(obj$coefficients)
  VC <- vcov(obj)
    
  # Draw values from the multivariate normal distribution
  Drawn <- rmultnorm(n = nsim, mu = Coef, vmat = VC)
  DrawnDF <- data.frame(Drawn)
  dfn <- names(DrawnDF)
 
  # If all values aren't set for calculating the hazard rate
  if (!isTRUE(means)){

    # Extract simulations for variables of interest
    bpos <- match(b, dfn)
    btvcpos <- match(btvc, dfn)
    
    Drawn <- data.frame(Drawn[, c(bpos, btvcpos)])
    Drawn$ID <- 1:nsim

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
        Xs <- data.frame(Xj)
        names(Xs) <- c("Xj")
        Xs$Comparison <- paste(Xs[, 1])
        Simb <- merge(TVSim, Xs)
        Simb$QI <- exp(Simb$CombCoef * Simb$Xj)
    } else if (qi == "First Difference"){
      if (length(Xj) != length(Xl)){
        stop("Xj and Xl must be the same length.")
      } else {
        TVSim$QI <- exp(TVSim$CombCoef)
        Xs <- data.frame(Xj, Xl)
        Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
        Simb <- merge(TVSim, Xs)
        Simb$QI <- (exp((Simb$Xj - Simb$Xl) * Simb$CombCoef) - 1) * 100
      }
    } else if (qi == "Hazard Ratio"){
       if (length(Xj) != length(Xl)){
        stop("Xj and Xl must be the same length.")
      } else {
        Xs <- data.frame(Xj, Xl)
        Xs$Comparison <- paste(Xs[, 1], "vs.", Xs[, 2])
        Simb <- merge(TVSim, Xs)
        Simb$QI <- exp((Simb$Xj - Simb$Xl) * Simb$CombCoef)
      }
    } else if (qi == "Hazard Rate"){
        Xl <- NULL
        message("Xl is ignored.")

        if (isTRUE(MeansMessage)){
          message("All variables values other than b are fitted at 0.")
        } 
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
        Simb$QI <- Simb$hazard * Simb$HR 
    }
  }

  # If the user wants to calculate Hazard Rates using means for fitting all covariates other than b.
  else if (isTRUE(means)){
    Xl <- NULL
    message("Xl ignored")

    # Set all values of b at means for data used in the analysis
    NotB <- setdiff(names(obj$means), b)
    MeanValues <- data.frame(obj$means)
    FittedMeans <- function(Z){
      ID <- 1:nsim
      Temp <- data.frame(ID)
      for (i in Z){
        BarValue <- MeanValues[i, ]
        DrawnCoef <- DrawnDF[, i]
        FittedCoef <- outer(DrawnCoef, BarValue)
        FCMolten <- data.frame(melt(FittedCoef))
        Temp <- cbind(Temp, FCMolten[,3])
      }
      Names <- c("ID", Z)
      names(Temp) <- Names
      Temp <- Temp[, -1]
      return(Temp)
    }
    FittedComb <- data.frame(FittedMeans(NotB)) 
    ExpandFC <- do.call(rbind, rep(list(FittedComb), length(Xj)))

    # Set fitted values for Xj
    bpos <- match(b, dfn)
    Simb <- data.frame(DrawnDF[, bpos])

    Xs <- data.frame(Xj)
    Xs$HRValue <- paste(Xs[, 1])

    Simb <- merge(Simb, Xs)
    Simb$CombB <- Simb[, 1] * Simb[, 2]
    Simb <- Simb[, 2:4]

    Simb <- cbind(Simb, ExpandFC)
    Simb$Sum <- rowSums(Simb[, c(-1, -2)])
    Simb$HR <- exp(Simb$Sum)
    Simb <- Simb[, c("HRValue", "HR", "Xj")]

    bfit <- basehaz(obj)
    bfit$FakeID <- 1
    Simb$FakeID <- 1
    bfitDT <- data.table(bfit, key = "FakeID", allow.cartesian = TRUE)
    SimbDT <- data.table(Simb, key = "FakeID", allow.cartesian = TRUE)
    SimbCombDT <- SimbDT[bfitDT, allow.cartesian = TRUE]
    Simb <- data.frame(SimbCombDT)
    Simb$QI <- Simb$hazard * Simb$HR 
  }

  # Drop simulations outside of 'confidence bounds'
  SubVar <- c("time", "Xj")
  
  SimbPerc <- IntervalConstrict(Simb = Simb, SubVar = SubVar, qi = qi,
                                QI = QI, spin = spin, ci = ci)

  # Create real time variable
  if (tfun == "linear"){
    SimbPerc$RealTime <- SimbPerc$tf
  } else if (tfun == "log"){
    SimbPerc$RealTime <- exp(SimbPerc$tf)
  } else if (tfun == "power"){
    SimbPerc$RealTime <- SimbPerc$tf^(1/pow)
  }
  
  # Final clean up
  class(SimbPerc) <- c("simtvc", qi)
  SimbPerc
}
