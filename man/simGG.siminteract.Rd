\name{simGG.siminteract}
\alias{simGG.siminteract}
\title{Plot simulated linear multiplicative interactions from Cox Proportional
Hazards Models}
\usage{
\method{simGG}{siminteract}(obj, from = NULL, to = NULL, xlab = NULL,
  ylab = NULL, title = NULL, method = "auto", spalette = "Set1",
  legend = "legend", leg.name = "", lcolour = "#2B8CBE", lsize = 1,
  pcolour = "#A6CEE3", psize = 1, alpha = 0.2, type = "points", ...)
}
\arguments{
  \item{obj}{a \code{siminteract} class object}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

  \item{from}{numeric time to start the plot from. Only
  relevant if \code{qi = "Hazard Rate"}.}

  \item{to}{numeric time to plot to. Only relevant if
  \code{qi = "Hazard Rate"}.}

  \item{title}{the plot's main title.}

  \item{method}{what type of smoothing method to use to
  summarize the center of the simulation distribution.}

  \item{spalette}{colour palette for when there are
  multiple sets of comparisons to plot. Not relevant for
  \code{qi = "Marginal Effect"}. Default palette is
  \code{"Set1"}. See \code{\link{scale_colour_brewer}}.}

  \item{legend}{specifies what type of legend to include
  (if applicable). The default is \code{legend = "legend"}.
  To hide the legend use \code{legend = FALSE}. See the
  \code{\link{discrete_scale}} for more details.}

  \item{leg.name}{name of the legend (if applicable).}

  \item{lcolour}{character string colour of the smoothing
  line. The default is hexadecimal colour \code{lcolour =
  '#2B8CBE'}. Only relevant if \code{qi = "Marginal
  Effect"}.}

  \item{lsize}{size of the smoothing line. Default is 1.
  See \code{\link{ggplot2}}.}

  \item{pcolour}{character string colour of the simulated
  points or ribbons (when there are not multiple sets of
  simulations). Default is hexadecimal colour \code{pcolour
  = '#A6CEE3'}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{alpha}{numeric. Alpha (e.g. transparency) for the
  points, lines, or ribbons. Default is \code{alpha = 0.2}.
  See \code{\link{ggplot2}}. Note, if \code{type = "lines"}
  or \code{type = "points"} then \code{alpah} sets the
  maximum value per line or point at the center of the
  distribution. Lines or points further from the center are
  more transparent the further they get from the middle.}

  \item{type}{character string. Specifies how to plot the
  simulations. Can be \code{points}, \code{lines}, or
  \code{ribbons}. If points then each simulation value will
  be plotted. If \code{lines} is chosen then each
  simulation is plotted using a different line. Note: any
  simulation with a value along its length that is outside
  of the specified central interval will be dropped. This
  is to create a smooth plot. If \code{type = "ribbons"} a
  plot will be created with shaded areas ('ribbons') for
  the minimum and maximum simulation values (i.e. the
  middle interval set with \code{qi} in
  \code{\link{coxsimSpline}}) as well as the central 50
  percent of this area. It also plots a line for the median
  value of the full area, so values in \code{method} are
  ignored. One of the key advantages of using ribbons
  rather than points is that it creates plots with smaller
  file sizes.}

  \item{...}{Additional arguments. (Currently ignored.)}
}
\value{
a \code{gg} \code{ggplot} class object
}
\description{
\code{simGG.siminteract} uses \link{ggplot2} to plot the
quantities of interest from \code{siminteract} objects,
including marginal effects, first differences, hazard
ratios, and hazard rates.
}
\details{
Uses \link{ggplot2} to plot the quantities of interest from
\code{siminteract} objects, including marginal effects,
first differences, hazard ratios, and hazard rates. If
there are multiple strata, the quantities of interest will
be plotted in a grid by strata. Note: A dotted line is
created at y = 1 (0 for first difference), i.e. no effect,
for time-varying hazard ratio graphs. No line is created
for hazard rates.

Note: if \code{qi = "Hazard Ratio"} or \code{qi = "First
Difference"} then you need to have choosen more than one
fitted value for \code{X1} in \code{\link{coxsimInteract}}.
}
\examples{
# Load Carpenter (2002) data
data("CarpenterFdaData")

# Load survival package
library(survival)

# Run basic model
M1 <- coxph(Surv(acttime, censor) ~ lethal*prevgenx,
            data = CarpenterFdaData)

# Simulate Marginal Effect of lethal for multiple values of prevgenx
Sim1 <- coxsimInteract(M1, b1 = "lethal", b2 = "prevgenx",
                     X2 = seq(2, 115, by = 2), nsim = 100)

# Plot quantities of interest
simGG(Sim1)

\dontrun{
# Change the order of the covariates to make a more easily
# interpretable hazard ratio graph.
M2 <- coxph(Surv(acttime, censor) ~ prevgenx*lethal,
             data = CarpenterFdaData)

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

# Simulate Hazard Rate
Sim4 <- coxsimInteract(M2, b1 = "prevgenx", b2 = "lethal",
                       X1 = 100, X2 = 1, qi = "Hazard Rate")

# Plot quantities of interest
simGG(Sim1, xlab = "\\nprevgenx",
       ylab = "Marginal Effect of lethal\\n")
simGG(Sim2, ribbons = TRUE)
simGG(Sim3)
simGG(Sim4, to = 150, ribbons = TRUE, legend = FALSE)
}
}
\references{
Brambor, Thomas, William Roberts Clark, and Matt Golder.
2006. ''Understanding Interaction Models: Improving
Empirical Analyses.'' Political Analysis 14(1): 63-82.

Keele, Luke. 2010. ''Proportionally Difficult: Testing for
Nonproportional Hazards in Cox Models.'' Political Analysis
18(2): 189-205.

Carpenter, Daniel P. 2002. ''Groups, the Media, Agency
Waiting Costs, and FDA Drug Approval.'' American Journal of
Political Science 46(3): 490-505.
}
\seealso{
\code{\link{coxsimInteract}},
\code{\link{simGG.simlinear}}, and \code{\link{ggplot2}}
}

