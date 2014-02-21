\name{tvc}
\alias{tvc}
\title{Create a time interaction variable}
\usage{
tvc(data, b, tvar, tfun = "linear", pow = NULL)
}
\arguments{
  \item{data}{a data frame}

  \item{b}{the non-time interacted variable's name}

  \item{tvar}{the time variable's name}

  \item{tfun}{function of time that btvc was multiplied by.
  Default is \code{tfun = "linear"}. Can also be \code{tfun
  = "log"} (natural log) and \code{tfun = "power"}. If
  \code{tfun = "power"} then the pow argument needs to be
  specified also.}

  \item{pow}{if \code{tfun = "power"}, then use pow to
  specify what power to raise the time interaction to.}
}
\value{
a vector
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
data("GolubEUPData")

# Expand data into equally spaced time intervals
GolubEUPDataExp <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                     Time = 'begin', Time2 = 'end', event = 'event')

# Create natural log time interaction with the qmv variable
GolubEUPData$Lqmv <- tvc(GolubEUPData, b = "qmv",
                       tvar = "end", tfun = "log")
}
\seealso{
\code{\link{SurvExpand}}, \code{\link{simGG.simtvc}},
\code{coxsimtvc}, \code{\link{survival}}, and
\code{\link{coxph}}
}
\keyword{utilities}

