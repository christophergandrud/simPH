\name{SurvExpand}
\alias{SurvExpand}
\title{Convert a data frame of non-equal interval continuous observations into
equal interval continuous observations}
\usage{
SurvExpand(data, GroupVar, Time, Time2, event, messages = TRUE)
}
\arguments{
  \item{data}{a data frame.}

  \item{GroupVar}{a character string naming the unit
  grouping variable.}

  \item{Time}{a character string naming the variable with
  the interval start time.}

  \item{Time2}{a character string naming the variable with
  the interval end time.}

  \item{event}{a character string naming the event
  variable. Note: must be numeric with 0 indicating no
  event.}

  \item{messages}{logical indicating if you want messages
  returned while the function is working.}
}
\value{
Returns a data frame where observations have been expanded
into equally spaced time intervals.
}
\description{
\code{SurvExpand} convert a data frame of non-equal
interval continuous observations into equal interval
continuous observations. This is useful for creating
time-interactions with \code{\link{tvc}}.
}
\details{
The function primarily prepares data from the creation of
accurate time-interactions with the \code{\link{tvc}}
command. Note: the function will work best if your original
time intervals are recorded in whole numbers. It also
currently does not support repeated events data.
}
\examples{
\dontrun{
# Load Golub & Steunenberg (2007) Data
data("GolubEUPData")

# Subset PURELY TO SPEED UP THE EXAMPLE
GolubEUPData <- GolubEUPData[1:500, ]

# Expand data
GolubEUPDataExpand <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                       Time = 'begin', Time2 = 'end', event = 'event')
}
}
\seealso{
\code{\link{tvc}}
}
\keyword{utilities}

