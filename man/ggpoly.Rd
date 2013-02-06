\name{ggpoly}
\alias{ggpoly}
\title{Plot simulated polynomial hazards.}
\usage{
  ggpoly(obj, xlab = NULL, ylab = NULL, title = NULL,
    xbreaks = NULL, xlabels = NULL, smoother = "auto",
    colour = "#A6CEE3", spalette = "Set1", leg.name = "",
    lsize = 2, psize = 1, palpha = 0.1, ...)
}
\arguments{
  \item{obj}{a simpoly object}
}
\description{
  \code{ggpoly} uses ggplot2 to plot simulated relative
  hazards from a simpoly class object.
}

