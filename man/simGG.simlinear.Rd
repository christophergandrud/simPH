\name{simGG.simlinear}
\alias{simGG.simlinear}
\title{Plot simulated linear time-constant quantities of interest from Cox Proportional Hazards Models}
\usage{
  \method{simGG}{simlinear} (obj, from = NULL, to = NULL,
    xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", spalette = "Set1", leg.name = "",
    lcolour = "#2B8CBE", lsize = 1, pcolour = "#A6CEE3",
    psize = 1, alpha = 0.1, ribbons = FALSE, ...)
}
\arguments{
  \item{obj}{a \code{simlinear} class object.}

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

  \item{leg.name}{name of the stratified hazard rates
  legend. Only relevant if \code{qi = "Hazard Rate"}.}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "First
  Difference"}.}

  \item{lsize}{size of the smoothing line. Default is 1.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour \code{pcolour = '#A6CEE3'}. Only relevant if
  \code{qi = "First Difference"}.}

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
  \code{\link{coxsimLinear}}) as well as the central 50
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
  \code{simGG.simlinear} uses \link{ggplot2} to plot the
  quantities of interest from \code{simlinear} objects,
  including relative hazards, first differences, hazard
  ratios, and hazard rates.
}
\details{
  Uses \link{ggplot2} to plot the quantities of interest
  from \code{simlinear} objects, including relative
  hazards, first differences, hazard ratios, and hazard
  rates. If there are multiple strata, the quantities of
  interest will be plotted in a grid by strata. Note: A
  dotted line is created at y = 1 (0 for first difference),
  i.e. no effect, for time-varying hazard ratio graphs. No
  line is created for hazard rates.
}
\examples{
## dontrun
# Load survival package
# library(survival)
# Load Carpenter (2002) data
# data("CarpenterFdaData")

# Estimate survival model
# M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
#            deathrt1 + acutediz + hosp01  + hhosleng +
#            mandiz01 + femdiz01 + peddiz01 + orphdum +
#            vandavg3 + wpnoavg3 + condavg3 + orderent +
#            stafcder, data = CarpenterFdaData)

# Simulate and plot Hazard Ratios for stafcder variable
# Sim1 <- coxsimLinear(M1, b = "stafcder",
#						qi = "Hazard Ratio",
#						Xj = seq(1237, 1600, by = 2), spin = TRUE)
# simGG(Sim1)

## dontrun
# Simulate and plot Hazard Rate for stafcder variable
# Sim2 <- coxsimLinear(M1, b = "stafcder", nsim = 100,
#						qi = "Hazard Rate", Xj = c(1237, 1600))
# simGG(Sim2, ribbons = TRUE)
}
\references{
  Licht, Amanda A. 2011. ''Change Comes with Time:
  Substantive Interpretation of Nonproportional Hazards in
  Event History Analysis.'' Political Analysis 19: 227-43.

  Keele, Luke. 2010. ''Proportionally Difficult: Testing
  for Nonproportional Hazards in Cox Models.'' Political
  Analysis 18(2): 189-205.

  Carpenter, Daniel P. 2002. ''Groups, the Media, Agency
  Waiting Costs, and FDA Drug Approval.'' American Journal
  of Political Science 46(3): 490-505.
}
\seealso{
  \code{\link{coxsimLinear}}, \code{\link{simGG.simtvc}},
  and \code{\link{ggplot2}}
}

