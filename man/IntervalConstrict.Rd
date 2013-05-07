\name{IntervalConstrict}
\alias{IntervalConstrict}
\title{Constrict simulations to a defined interval.}
\usage{
  IntervalConstrict(Simb = Simb, SubVar = SubVar, qi = qi,
    QI = QI, spin = FALSE, ci = 0.95)
}
\arguments{
  \item{Simb}{character string naming the data frame with
  the simulations.}

  \item{SubVar}{character vector the variable names to
  subset the simulations by.}

  \item{qi}{character vector naming the type of quantity of
  interest.}

  \item{QI}{character string labeling the quantitiy of
  interest values.}

  \item{spin}{logical for whether or not to use the
  shortest probability interval or the central interval.}

  \item{ci}{numeric confidence interval measure.}
}
\description{
  \code{IntervalConstrict} is an internal function to
  constrict a set of simulations to user defined interval.
}

