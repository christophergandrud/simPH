\name{simGG.simtvc}
\alias{simGG.simtvc}
\title{Plot simulated time-varying hazard ratios or stratified time-varying hazard rates from a simtvc class object using ggplot2}
\usage{
  \method{simGG}{simtvc} (obj, from = NULL, to = NULL,
    xlab = NULL, ylab = NULL, title = NULL,
    smoother = "auto", spalette = "Set1", leg.name = "",
    lcolour = "#2B8CBE", lsize = 2, pcolour = "#A6CEE3",
    psize = 1, palpha = 0.1)
}
\arguments{
  \item{obj}{a simtvc class object}

  \item{from}{numeric time to start the plot from.}

  \item{to}{numeric time to plot to.}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

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
}
\value{
  a ggplot2 object
}
\description{
  \code{simGG.simtvc} uses ggplot2 to plot the simulated
  hazards from a \code{simtvc} class object using ggplot2.
}
\details{
  Plots either a time varying hazard ratio or the hazard
  rates for multiple strata. Currently the strata legend
  labels need to be changed manually (see
  \code{\link{revalue}} in the \link{plyr} package) in the
  \code{simtvc} object with the \code{strata} component.
  Also, currently the x-axis tick marks and break labels
  must be adjusted manually for non-linear functions of
  time. Note: A dotted line is created at y = 1 (0 for
  first difference), i.e. no effect, for time-varying
  hazard ratio graphs. No line is created for hazard rates.
}
\examples{
# Load Golub & Steunenberg (2007) Data
data("GolubEUPData")

# Load survival package
library(survival)

# Create natural log time interactions
Golubtvc <- function(x){
  assign(paste0("l", x), tvc(GolubEUPData, b = x, tvar = "end", tfun = "log"))
}

GolubEUPData$Lcoop <-Golubtvc("coop")
GolubEUPData$Lqmv <- Golubtvc("qmv")
GolubEUPData$Lbacklog <- Golubtvc("backlog")
GolubEUPData$Lcodec <- Golubtvc("codec")
GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
GolubEUPData$Lthatcher <- Golubtvc("thatcher")

# Run Cox PH Model
M1 <- coxph(Surv(begin, end, event) ~
            qmv + qmvpostsea + qmvpostteu +
            coop + codec + eu9 + eu10 + eu12 +
            eu15 + thatcher + agenda + backlog +
            Lqmv + Lqmvpostsea + Lcoop + Lcodec +
            Lthatcher + Lbacklog,
         data = GolubEUPData,
         ties = "efron")

# Create simtvc object for Relative Hazard
Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
                   tfun = "log", from = 80, to = 2000,
                   Xj = 1, by = 15, ci = 0.99)

# Create simtvc object for First Difference
Sim2 <- coxsimtvc(obj = M1, b = "qmv", btvc = "Lqmv",
                 qi = "First Difference", Xj = 1,
                 tfun = "log", from = 80, to = 2000,
                 by = 15, ci = 0.95)

# Create simtvc object for Hazard Ratio
Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
                  qi = "Hazard Ratio", Xj = c(191, 229),
                  Xl = c(0, 0),
                  tfun = "log", from = 80, to = 2000,
                  by = 15, ci = 0.99)

# Create plots
simGG(Sim1)
simGG(Sim2)
simGG(Sim3, leg.name = "Comparision", from = 1200)
}
\references{
  Licht, Amanda A. 2011. ''Change Comes with Time:
  Substantive Interpretation of Nonproportional Hazards in
  Event History Analysis.'' Political Analysis 19: 227-43.
}

