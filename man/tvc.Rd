\name{tvc}
\alias{tvc}
\title{Create a time interaction variable}
\usage{
tvc(data, b, tvar, tfun = "linear", pow = NULL, vector = FALSE)
}
\arguments{
  \item{data}{a data frame}

  \item{b}{the non-time interacted variable's name. Either
  a single value or a vector of variable names can be
  entered.}

  \item{tvar}{the time variable's name}

  \item{tfun}{function of time that btvc was multiplied by.
  Default is \code{tfun = "linear"}. Can also be \code{tfun
  = 'log'} (natural log) and \code{tfun = 'power'}. If
  \code{tfun = 'power'} then the pow argument needs to be
  specified also.}

  \item{pow}{if \code{tfun = 'power'}, then use pow to
  specify what power to raise the time interaction to.}

  \item{vector}{logical. Whether or not to return one
  vector a or a data frame. Can only be used if only one
  \code{b} is included.}
}
\value{
a data frame or vector. If a data frame is returned it will
include all of the original variables as well as the
interactions denoted by a variable name
'\code{bn}_\code{tfun}', where \code{bn} is one variable
name from \code{b} and \code{tfun} as entered into the
function.
}
\description{
\code{tvc} creates a time interaction variable that can be
used in a coxph model (or any other model with time
interactions)
}
\details{
Interacts a variable with a specified function of time.
Possible functions of time include \code{'linear'}, natural
\code{'log'}, and exponentiated (\code{'power'}).
}
\examples{
# Load Golub & Steunenberg (2007) Data
data('GolubEUPData')

# Subset PURELY TO SPEED UP THE EXAMPLE
GolubEUPData <- GolubEUPData[1:500, ]

# Expand data into equally spaced time intervals
GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                  Time = 'begin', Time2 = 'end', event = 'event')

# Create natural log time interaction with the qmv variable
GolubEUPData$Lqmv <- tvc(GolubEUPData, b = 'qmv', tvar = 'end', tfun = 'log',
                         vector = TRUE)

# Create interactions for a vector of variables
BaseVars <- c('qmv', 'backlog', 'coop', 'codec', 'qmvpostsea', 'thatcher')

Test <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')
}
\seealso{
\code{\link{SurvExpand}}, \code{\link{simGG.simtvc}},
\code{\link{coxsimtvc}}, \code{\link{survival}}, and
\code{\link{coxph}}
}
\keyword{utilities}

