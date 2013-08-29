\name{simGG.simspline}
\alias{simGG.simspline}
\title{Plot simulated penalised spline hazards from Cox Proportional Hazards Models}
\usage{
  \method{simGG}{simspline} (obj, FacetTime = NULL,
    from = NULL, to = NULL, xlab = NULL, ylab = NULL,
    zlab = NULL, title = NULL, smoother = "auto",
    leg.name = "", lcolour = "#2B8CBE", lsize = 2,
    pcolour = "#A6CEE3", psize = 1, palpha = 0.1,
    surface = TRUE, fit = "linear", ribbons = FALSE, ...)
}
\arguments{
  \item{obj}{a \code{simspline} class object}

  \item{FacetTime}{a numeric vector of points in time where
  you would like to plot Hazard Rates in a facet grid. Only
  relevant if \code{qi == 'Hazard Rate'}. Note: the values
  of Facet Time must exactly match values of the
  \code{time} element of \code{obj}.}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

  \item{zlab}{a label for the plot's z-axis. Only relevant
  if \code{qi = "Hazard Rate"} and \code{FacetTime ==
  NULL}.}

  \item{from}{numeric time to start the plot from. Only
  relevant if \code{qi = "Hazard Rate"}.}

  \item{to}{numeric time to plot to. Only relevant if
  \code{qi = "Hazard Rate"}.}

  \item{title}{the plot's main title.}

  \item{smoother}{what type of smoothing line to use to
  summarize the plotted coefficient.}

  \item{leg.name}{name of the stratified hazard rates
  legend. Only relevant if \code{qi = "Hazard Rate"}.}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "Relative
  Hazard"} or \code{qi = "First Difference"}.}

  \item{lsize}{size of the smoothing line. Default is 2.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour \code{pcolour = '#A6CEE3'}. Only relevant if
  \code{qi = "Relative Hazard"} or \code{qi = "First
  Difference"}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{palpha}{point alpha (e.g. transparency). Default is
  \code{palpha = 0.05}. See \code{\link{ggplot2}}.}

  \item{surface}{plot surface. Default is \code{surface =
  TRUE}. Only relevant if \code{qi == 'Relative Hazard'}
  and \code{FacetTime = NULL}.}

  \item{fit}{one or more of \code{"linear"},
  \code{"quadratic"}, \code{"smooth"}, \code{"additive"};
  to display fitted surface(s); partial matching is
  supported e.g., \code{c("lin", "quad")}. Only relevant if
  \code{qi == 'Relative Hazard'} and \code{FacetTime =
  NULL}.}

  \item{ribbons}{logical specifies whether or not to use
  summary ribbons of the simulations rather than plotting
  every simulation value as a point. If \code{lines = TRUE}
  a plot will be created with shaded areas ('ribbons') for
  the minimum and maximum simulation values (i.e. the
  middle interval set with \code{qi} in
  \code{\link{coxsimSpline}}) as well as the central 50% of
  this area. It also plots a line for the median value of
  the full area, so values in \code{smoother} are ignored.
  One of the key advantages of using ribbons rather than
  points is that it creates plots with smaller file sizes.}

  \item{...}{Additional arguments. (Currently ignored.)}
}
\value{
  a \code{gg} \code{ggplot} class object. See
  \code{\link{scatter3d}} for values from \code{scatter3d}
  calls.
}
\description{
  \code{simGG.simspline} uses \link{ggplot2} and
  \link{scatter3d} to plot quantities of interest from
  \code{simspline} objects, including relative hazards,
  first differences, hazard ratios, and hazard rates.
}
\details{
  Uses \code{ggplot2} and \code{scatter3d} to plot the
  quantities of interest from \code{simspline} objects,
  including relative hazards, first differences, hazard
  ratios, and hazard rates. If currently does not support
  hazard rates for multiple strata.

  It can graph hazard rates as a 3D plot using
  \code{\link{scatter3d}} with the dimensions: Time, Hazard
  Rate, and the value of \code{Xj}. Ribbon plots are not
  available with 3D plots. You can also choose to plot
  hazard rates for a range of values of \code{Xj} in two
  dimensional plots at specific points in time. Each plot
  is arranged in a facet grid.

  Note: A dotted line is created at y = 1 (0 for first
  difference), i.e. no effect, for time-varying hazard
  ratio graphs. No line is created for hazard rates.
}
\examples{
## dontrun
# Load Carpenter (2002) data
# data("CarpenterFdaData")

# Load survival package
# library(survival)

# Run basic model
# From Keele (2010) replication data
# M1 <- coxph(Surv(acttime, censor) ~  prevgenx + lethal + deathrt1 +
#				acutediz + hosp01  + pspline(hospdisc, df = 4) +
#				pspline(hhosleng, df = 4) + mandiz01 + femdiz01 + peddiz01 +
#				orphdum + natreg + vandavg3 + wpnoavg3 +
#				pspline(condavg3, df = 4) + pspline(orderent, df = 4) +
#				pspline(stafcder, df = 4), data = CarpenterFdaData)

# Simulate Relative Hazards for orderent
# Sim1 <- coxsimSpline(M1, bspline = "pspline(stafcder, df = 4)",
#                    bdata = CarpenterFdaData$stafcder,
#                    qi = "Hazard Ratio",
#                    Xj = seq(1100, 1700, by = 10),
#                    Xl = seq(1099, 1699, by = 10), spin = TRUE)

# Simulate Hazard Rate for orderent
# Sim2 <- coxsimSpline(M1, bspline = "pspline(orderent, df = 4)",
#                    bdata = CarpenterFdaData$orderent,
#                    qi = "Hazard Rate",
#                    Xj = seq(1, 30, by = 2), ci = 0.9, nsim = 10)

# Plot relative hazard
# simGG(Sim1, palpha = 1)

# 3D plot hazard rate
# simGG(Sim2, zlab = "orderent", fit = "quadratic")

# Create a time grid plot
# Find all points in time where baseline hazard was found
# unique(Sim2$Time)

# Round time values so they can be exactly matched with FacetTime
# Sim2$Time <- round(Sim2$Time, digits = 2)

# Create plot
# simGG(Sim2, FacetTime = c(6.21, 25.68, 100.64, 202.36), ribbons = TRUE, palpha = 0.5)

# Simulated Fitted Values of stafcder
# Sim3 <- coxsimSpline(M1, bspline = "pspline(stafcder, df = 4)",
#                    bdata = CarpenterFdaData$stafcder,
#                    qi = "Hazard Ratio",
#                    Xj = seq(1100, 1700, by = 10),
#                    Xl = seq(1099, 1699, by = 10), ci = 0.90)

# Plot simulated Hazard Ratios
# simGG(Sim3, xlab = "\\nFDA Drug Review Staff", palpha = 0.2)
}
\seealso{
  \code{\link{coxsimLinear}}, \code{\link{simGG.simtvc}},
  \code{\link{ggplot2}}, and \code{\link{scatter3d}}
}

