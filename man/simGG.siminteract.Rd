\name{simGG.siminteract}
\alias{simGG.siminteract}
\title{Plot simulated linear multiplicative interactions.}
\usage{
  \method{simGG}{siminteract} (obj, from = NULL, to = NULL,
    xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", spalette = "Set1", leg.name = "",
    lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3",
    psize = 1, palpha = 0.1)
}
\arguments{
  \item{obj}{a simlinear object}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

  \item{from}{numeric time to start the plot from. Only
  relevant if \code{qi = "Hazard Rate"}.}

  \item{to}{numeric time to plot to. Only relevant if
  \code{qi = "Hazard Rate"}.}

  \item{title}{the plot's main title}

  \item{smoother}{what type of smoothing line to use to
  summarize the plotted coefficient}

  \item{spalette}{colour palette. Not relevant for \code{qi
  = "Marginal Effect"}. Default palette is \code{"Set1"}.
  See \code{\link{scale_colour_brewer}}.}

  \item{leg.name}{name of the stratified hazard rates
  legend. Only relevant if \code{qi = "Hazard Rate"}.}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "Marginal
  Effect"}.}

  \item{lsize}{size of the smoothing line. Default is 2.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour \code{pcolour = '#A6CEE3'}. Only relevant if
  \code{qi = "Marginal Effect"}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{palpha}{point alpha (e.g. transparency). Default is
  \code{palpha = 0.05}. See \code{\link{ggplot2}}.}
}
\value{
  a ggplot2 object
}
\description{
  Uses ggplot2 to plot the quantities of interest from
  \code{siminteract} objects, including marginal effects,
  first differences, hazard ratios, and hazard rates. If
  there are multiple strata, the quantities of interest
  will be plotted in a grid by strata. Note: A dotted line
  is created at y = 1 (0 for first difference), i.e. no
  effect, for time-varying hazard ratio graphs. No line is
  created for hazard rates.

  Note: if \code{qi = "Hazard Ratio"} or \code{qi = "First
  Difference"} then you need to have choosen more than one
  fitted value for \code{X1} in
  \code{\link{coxsimInteract}}.
}
\details{
  \code{simGG.siminteract} uses ggplot2 to plot the
  quantities of interest from \code{siminteract} objects,
  including marginal effects, first differences, hazard
  ratios, and hazard rates.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx, data = CarpenterFdaData)

# Simulate Marginal Effect of lethal for multiple values of prevgenx
Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx", X2 = seq(2, 115, by = 2))
# Change the order of the covariates to make a more easily
# interpretable hazard ratio graph.
M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal, data = CarpenterFdaData)

# Simulate Hazard Ratio of lethal for multiple values of prevgenx
Sim2 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
                    X1 = seq(2, 115, by = 2),
                    X2 = c(0, 1),
                    qi = "Hazard Ratio", ci = 0.9)

# Simulate First Difference
Sim3 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
                       X1 = seq(2, 115, by = 2),
                       X2 = c(0, 1),
                       qi = "First Difference", spin = TRUE)

# Plot quantities of interest
simGG(Sim1, xlab = "\\nprevgenx", ylab = "Marginal Effect of lethal\\n")
simGG(Sim2)
simGG(Sim3)
}
\references{
  Brambor, Thomas, William Roberts Clark, and Matt Golder.
  2006. ''Understanding Interaction Models: Improving
  Empirical Analyses.'' Political Analysis 14(1): 63-82.

  Keele, Luke. 2010. ''Proportionally Difficult: Testing
  for Nonproportional Hazards in Cox Models.'' Political
  Analysis 18(2): 189-205.

  Carpenter, Daniel P. 2002. ''Groups, the Media, Agency
  Waiting Costs, and FDA Drug Approval.'' American Journal
  of Political Science 46(3): 490-505.
}
\seealso{
  \code{\link{coxsimInteract}},
  \code{\link{simGG.simlinear}}, and \code{\link{ggplot2}}
}

