\name{simGG.simpoly}
\alias{simGG.simpoly}
\title{Plot simulated polynomial quantities of interest from Cox Proportional Hazards Models}
\usage{
  \method{simGG}{simpoly} (obj, from = NULL, to = NULL,
    xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", spalette = "Set1",
    legend = "legend", leg.name = "", lcolour = "#2B8CBE",
    lsize = 1, pcolour = "#A6CEE3", psize = 1, alpha = 0.1,
    ribbons = FALSE, ...)
}
\arguments{
  \item{obj}{a \code{simpoly} class object.}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

  \item{from}{numeric time to start the plot from. Only
  relevant if \code{qi = "Hazard Rate"}.}

  \item{to}{numeric time to plot to. Only relevant if
  \code{qi = "Hazard Rate"}.}

  \item{title}{the plot's main title.}

  \item{smoother}{what type of smoothing line to use to
  summarize the plotted coefficient.}

  \item{spalette}{colour palette for use in \code{qi =
  "Hazard Rate"}. Default palette is \code{"Set1"}. See
  \code{\link{scale_colour_brewer}}.}

  \item{legend}{specifies what type of legend to include
  (if applicable). The default is \code{legend = "legend"}.
  To hide the legend use \code{legend = FALSE}. See the
  \code{\link{discrete_scale}} for more details.}

  \item{leg.name}{name of the legend (if applicable).}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "First
  Difference"}.}

  \item{lsize}{size of the smoothing line. Default is 1.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points or ribbons (when there are not multiple sets of
  simulations). Default is hexadecimal colour \code{pcolour
  = '#A6CEE3'}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{alpha}{point alpha (e.g. transparency) for the
  points or ribbons. Default is \code{alpha = 0.1}. See
  \code{\link{ggplot2}}.}

  \item{ribbons}{logical specifies whether or not to use
  summary ribbons of the simulations rather than plotting
  every simulation value as a point. If \code{ribbons =
  TRUE} a plot will be created with shaded areas
  ('ribbons') for the minimum and maximum simulation values
  (i.e. the middle interval set with \code{qi} in
  \code{\link{coxsimPoly}}) as well as the central 50
  percent of this area. It also plots a line for the median
  value of the full area, so values in \code{smoother} are
  ignored. One of the key advantages of using ribbons
  rather than points is that it creates plots with smaller
  file sizes.}

  \item{...}{Additional arguments. (Currently ignored.)}
}
\value{
  a \code{gg} \code{ggplot} class object
}
\description{
  \code{simGG.simpoly} uses \link{ggplot2} to plot
  simulated relative quantities of interest from a
  \code{simpoly} class object.
}
\details{
  Uses \link{ggplot2} to plot the quantities of interest
  from \code{simpoly} objects.
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
           pow = 3, Xj = seq(1, 150, by = 5), nsim = 100)

# dontrun
# Simulate simpoly Hazard Ratio with spin probibility interval
# Sim2 <- coxsimPoly(M1, b = "natreg", qi = "Hazard Ratio",
#           pow = 3, Xj = seq(1, 150, by = 5), spin = TRUE)

# Plot simulations
simGG(Sim1)
# dontrun
# simGG(Sim2)
}
\seealso{
  \code{\link{coxsimPoly}} and \code{\link{ggplot2}}
}

