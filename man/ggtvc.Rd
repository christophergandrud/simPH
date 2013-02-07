\name{ggtvc}
\alias{ggtvc}
\title{Plot simulated time-varying hazard ratios or stratified time-varying hazard rates from a simtvc class object using ggplot2}
\usage{
  ggtvc(obj, qi = "Relative Hazard", strata = FALSE,
    from = NULL, to = NULL, xlab = NULL, ylab = NULL,
    title = NULL, smoother = "auto", colour = "#A6CEE3",
    spalette = "Set1", leg.name = "", lsize = 2, psize = 1,
    palpha = 0.1, ...)
}
\arguments{
  \item{obj}{a simtvc class object}

  \item{qi}{character string indicating what quantity of
  interest you would like to calculate. Can be
  \code{'Relative Hazard'}, \code{'First Difference'}, or
  \code{'Hazard Ratio'}. Default is \code{qi = 'Relative
  Hazard'}.}

  \item{strata}{logical whether or not you would like to
  plot the hazard rate for the separate strata}

  \item{from}{numeric time to start the plot from.}

  \item{to}{numeric time to plot to.}

  \item{xlab}{a label for the plot's x-axis.}

  \item{ylab}{a label of the plot's y-axis. The default
  uses the value of \code{qi}.}

  \item{title}{the plot's main title}

  \item{smoother}{what type of smoothing line to use to
  summarize the plotted coefficient}

  \item{colour}{character string colour of the simulated
  points for relative hazards. Default is hexadecimal
  colour A6CEE3. Works if \code{strata = FALSE}.}

  \item{spalette}{colour palette for stratified hazard
  rates. Only works if \code{strata = TRUE}. Default
  palette is \code{"Set1"}. See
  \code{\link{scale_colour_brewer}}.}

  \item{leg.name}{name of the stratified hazard rates
  legend. Only works if \code{strata = TRUE}.}

  \item{lsize}{size of the smoothing line. Default is 2.
  See \code{\link{ggplot2}}.}

  \item{psize}{size of the plotted simulation points.
  Default is \code{psize = 1}. See \code{\link{ggplot2}}.}

  \item{palpha}{point alpha (e.g. transparency). Default is
  \code{palpha = 0.05}. See \code{\link{ggplot2}}.}

  \item{...}{other arguments passed to specific methods}
}
\value{
  a ggplot2 object
}
\description{
  \code{ggtvc} uses ggplot2 to plot the simulated hazards
  from a simtvc class object using ggplot2. Note: A dotted
  line is created at y = 1 (0 for first difference), i.e.
  no effect, for time-varying hazard ratio graphs.
}
\details{
  Plots either a time varying hazard ratio or the hazard
  rates for multiple strata. Currently the strata legend
  labels need to be changed manually (see
  \code{\link{revalue}} in the \link{plyr} package) in the
  \code{simtvc} object with the \code{strata} component.
  Also, currently the x-axis tick marks and break labels
  must be adjusted manually for non-linear functions of
  time.
}
\examples{
# Load Golub & Steunenberg (2007) Data
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
                   by = 15, ci = "99")

# Create simtvc object for First Difference
Sim2 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
                  qi = "First Difference",
                  tfun = "log", from = 80, to = 2000,
                  by = 15, ci = "99")

# Create simtvc object for Hazard Ratio
Sim3 <- coxsimtvc(obj = M1, b = "backlog", btvc = "Lbacklog",
                  qi = "Hazard Ratio", Xj = c(191, 229),
                  Xl = c(0, 0),
                  tfun = "log", from = 80, to = 2000,
                  by = 15, ci = "99")

# Create plots
ggtvc(Sim1, qi = "Relative Hazard")
ggtvc(Sim2, qi = "First Difference")
ggtvc(Sim3, qi = "Hazard Ratio", leg.name = "Comparision", from = 1200)
}
\references{
  Licht, Amanda A. 2011. “Change Comes with Time:
  Substantive Interpretation of Nonproportional Hazards in
  Event History Analysis.” Political Analysis 19: 227–43.
}

