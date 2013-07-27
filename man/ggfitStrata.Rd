\name{ggfitStrata}
\alias{ggfitStrata}
\title{Graph fitted stratified survival curves}
\usage{
  ggfitStrata(obj, byStrata = FALSE, xlab = "", ylab = "",
    title = "", lcolour = "#2C7FB8", rcolour = "#2C7FB8")
}
\arguments{
  \item{obj}{a \code{survfit} object.}

  \item{byStrata}{logical, whether or not you want to
  include all of the stratified survival curves on one plot
  or separate them into a grid arranged plot.}

  \item{xlab}{a label for the plot's x-axis}

  \item{ylab}{a label of the plot's y-axis}

  \item{title}{plot title.}

  \item{lcolour}{line color. Currently only works if
  \code{byStrata = TRUE}. The default it \code{lcolour =
  "#2C7FB8"} (a bluish hexadecimal colour)}

  \item{rcolour}{confidence bounds ribbon color. Either a
  single color or a vector of colours. The default it
  \code{lcolour = "#2C7FB8"} (a bluish hexadecimal colour)}
}
\description{
  This function largely improves
  \code{\link{plot.survfit}}. It plots the curves using
  \link{ggplot2} rather than base R graphics. One major
  advantage is the ability to split the survival curves
  into multiple plots and arrange them in a grid. This
  makes it easier to examine many strata at once. Otherwise
  they can be very bunched up.

  Note: the strata legend labels need to be changed
  manually (see \code{\link{revalue}}) in the
  \code{survfit} object with the \code{strata} component.
}
\details{
  \code{ggfitStrata} graphs fitted survival curves created
  with \code{\link{survfit}} using \link{ggplot2}.
}
\examples{
#dontrun
# Load survival
#library(survival)
# Subset data
#bladder1 <- bladder[bladder$enum < 5, ]
# Estimate coxph model
#M1 <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + cluster(id), bladder1)
# Survfit
#M1Fit <- survfit(M1)
# Plot strata in a grid
#ggfitStrata(M1Fit, byStrata = TRUE)
}
\seealso{
  \code{\link{survfit}}, \code{\link{ggplot2}} and
  \code{\link{strata}}
}

