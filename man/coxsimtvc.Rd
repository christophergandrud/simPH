\name{coxsimtvc}
\alias{coxsimtvc}
\title{Simulate time-interactive quantities of interest from Cox Proportional
Hazards models}
\usage{
coxsimtvc(obj, b, btvc, qi = "Relative Hazard", Xj = NULL, Xl = NULL,
  tfun = "linear", pow = NULL, nsim = 1000, from, to, by = 1,
  ci = 0.95, spin = FALSE)
}
\arguments{
  \item{obj}{a \code{\link{coxph}} fitted model object with
  a time interaction.}

  \item{b}{the non-time interacted variable's name.}

  \item{btvc}{the time interacted variable's name.}

  \item{qi}{character string indicating what quantity of
  interest you would like to calculate. Can be
  \code{'Relative Hazard'}, \code{'First Difference'},
  \code{'Hazard Ratio'}, \code{'Hazard Rate'}. Default is
  \code{qi = 'Relative Hazard'}. If \code{qi = 'First
  Difference'} or \code{qi = 'Hazard Ratio'} then you can
  set \code{Xj} and \code{Xl}.}

  \item{Xj}{numeric vector of fitted values for \code{b}.
  Must be the same length as \code{Xl} or \code{Xl} must be
  \code{NULL}.}

  \item{Xl}{numeric vector of fitted values for Xl. Must be
  the same length as \code{Xj}. Only applies if \code{qi =
  'First Difference'} or \code{qi = 'Hazard Ratio'}.}

  \item{nsim}{the number of simulations to run per point in
  time. Default is \code{nsim = 1000}.}

  \item{tfun}{function of time that btvc was multiplied by.
  Default is "linear". It can also be "log" (natural log)
  and "power". If \code{tfun = "power"} then the pow
  argument needs to be specified also.}

  \item{pow}{if \code{tfun = "power"}, then use pow to
  specify what power the time interaction was raised to.}

  \item{from}{point in time from when to begin simulating
  coefficient values}

  \item{to}{point in time to stop simulating coefficient
  values.}

  \item{by}{time intervals by which to simulate coefficient
  values.}

  \item{ci}{the proportion of simulations to keep. The
  default is \code{ci = 0.95}, i.e. keep the middle 95
  percent. If \code{spin = TRUE} then \code{ci} is the
  confidence level of the shortest probability interval.
  Any value from 0 through 1 may be used.}

  \item{spin}{logical, whether or not to keep only the
  shortest probability interval rather than the middle
  simulations. Currently not supported for hazard rates.}
}
\value{
a \code{simtvc} object
}
\description{
\code{coxsimtvc} simulates time-interactive relative
hazards, first differences, and hazard ratios from models
estimated with \code{\link{coxph}} using the multivariate
normal distribution. These can be plotted with
\code{\link{simGG}}.
}
\details{
Simulates time-varying relative hazards, first differences,
and hazard ratios using parameter estimates from
\code{coxph} models. Can also simulate hazard rates for
multiple strata.

Relative hazards are found using: \deqn{RH = e^{\beta_{1} +
\beta_{2}f(t)}}{RH = exp(\beta[1] + \beta[2] f(t))} where
\eqn{f(t)} is the function of time.

First differences are found using: \deqn{FD = (e^{(X_{j} -
X_{l}) (\beta_{1} + \beta_{2}f(t))} - 1) * 100}{FD =
(exp((X[j] - X[l])(\beta[1] + \beta[2]f(t)) - 1) * 100}
where \eqn{X_{j}}{X[j]} and \eqn{X_{l}}{X[l]} are some
values of \eqn{X} to contrast.

Hazard ratios are calculated using: \deqn{FD = e^{(X_{j} -
X_{l}) (\beta_{1} + \beta_{2}f(t))}}{FD = exp((X[j] -
X[l])(\beta[1] + \beta[2]f(t))} When simulating
non-stratifed time-varying harzards \code{coxsimtvc} uses
the point estimates for a given coefficient
\eqn{\hat{\beta}_{x}}{hat \beta[x]} and its time
interaction \eqn{\hat{\beta}_{xt}}{hat \beta[xt]} along
with the variance matrix (\eqn{\hat{V}(\hat{\beta})}{hat
V(hat \beta)}) estimated from a \code{coxph} model. These
are used to draw values of \eqn{\beta_{1}}{\beta[1]} and
\eqn{\beta_{2}}{\beta[2]} from the multivariate normal
distribution \eqn{N(\hat{\beta},\:
\hat{V}(\hat{\beta}))}{N(hat \beta, hat V (hat \beta))}.

When simulating stratified time-varying hazard rates
\eqn{H} for a given strata \eqn{k}, \code{coxsimtvc} uses:
\deqn{H_{kxt} = \hat{\beta}_{k0t}e^{\hat{\beta_{1}} +
\beta_{2}f(t)}}{H_{kxt} = hat \beta[k0t]exp(hat \beta[1] +
\beta[2]f(t))} The resulting simulation values can be
plotted using \code{\link{simGG}}.
}
\examples{
# Load Golub & Steunenberg (2007) Data
data("GolubEUPData")

# Load survival package
library(survival)

# Create natural log time interactions
Golubtvc <- function(x){
  assign(paste0("l", x), tvc(GolubEUPData, b = x,
         tvar = "end", tfun = "log"))
}
\dontrun{
# Expand data (not run to speed processing time, but should be run)
GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                     Time = 'begin', Time2 = 'end', event = 'event')
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
                   Xj = 1, by = 15, ci = 0.99, nsim = 100)

\dontrun{
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
                  by = 15, ci = 0.5)
}
}
\references{
Golub, Jonathan, and Bernard Steunenberg. 2007. ''How Time
Affects EU Decision-Making.'' European Union Politics 8(4):
555-66.

Licht, Amanda A. 2011. ''Change Comes with Time:
Substantive Interpretation of Nonproportional Hazards in
Event History Analysis.'' Political Analysis 19: 227-43.

King, Gary, Michael Tomz, and Jason Wittenberg. 2000.
''Making the Most of Statistical Analyses: Improving
Interpretation and Presentation.'' American Journal of
Political Science 44(2): 347-61.

Liu, Ying, Andrew Gelman, and Tian Zheng. 2013.
''Simulation-Efficient Shortest Probability Intervals.''
Arvix. \url{http://arxiv.org/pdf/1302.2142v1.pdf}.
}
\seealso{
\code{\link{simGG}}, \code{\link{survival}},
\code{\link{strata}}, and \code{\link{coxph}}
}

