\name{simGG.simpoly}
\alias{simGG.simpoly}
\title{Plot simulated polynomial quantities of interest from Cox Proportional Hazards Models.}
\usage{
  \method{simGG}{simpoly} (obj, from = NULL, to = NULL,
    xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", spalette = "Set1", leg.name = "",
    lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3",
    psize = 1, palpha = 0.1)
}
\arguments{
  \item{obj}{a simpoly object.}

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

  \item{spalette}{colour palette for use in \code{qi =
  "Hazard Rate"}. Default palette is \code{"Set1"}. See
  \code{\link{scale_colour_brewer}}.}

  \item{leg.name}{name of the stratified hazard rates
  legend. Only relevant if \code{qi = "Hazard Rate"}.}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "First
  Difference"}.}

  \item{lsize}{size of the smoothing line. Default is 2.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour \code{pcolour = '#A6CEE3'}. Only relevant if
  \code{qi = "First Difference"}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{palpha}{point alpha (e.g. transparency). Default is
  \code{palpha = 0.05}. See \code{\link{ggplot2}}.}
}
\value{
  a ggplot2 object.
}
\description{
  \code{simGG.simpoly} uses ggplot2 to plot simulated
  relative quantities of interest from a \code{simpoly}
  class object.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal + deathrt1 + acutediz +
       hosp01  + hhosleng + mandiz01 + femdiz01 + peddiz01 + orphdum +
       natreg + I(natreg^2) + I(natreg^3) + vandavg3 + wpnoavg3 +
       condavg3 + orderent + stafcder, data = CarpenterFdaData)

# Simulate simpoly First Difference
Sim1 <- coxsimPoly(M1, b = "natreg", qi = "First Difference",
           pow = 3, Xj = seq(1, 150, by = 5))

# Simulate simpoly Hazard Ratio with spin probibility interval
Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio",
           pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)

# Plot simulations
simGG(Sim1)
simGG(Sim2)
}
\seealso{
  \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
}

