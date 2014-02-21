#' Convert a data frame of non-equal interval continuous observations into equal interval continuous observations
#' 
#' \code{SurvExpand} convert a data frame of non-equal interval continuous observations into equal interval continuous observations. This is useful for creating time-interactions with \code{\link{tvc}}.
#' @param data a data frame.
#' @param GroupVar a character string naming the unit grouping variable.
#' @param Time a character string naming the variable with the interval start time.
#' @param Time2 a character string naming the variable with the interval end time.
#' @param event a character string naming the event variable. Note: must be numeric with 0 indicating no event.
#' @param messages logical indicating if you want messages returned while the function is working.
#' 
#' @return Returns a data frame where observations have been expanded into equally spaced time intervals.
#' 
#' @details The function primarily prepares data from the creation of accurate time-interactions with the \code{\link{tvc}} command.
#' Note: the function will work best if your original time intervals are recorded in whole numbers. It also currently does not support repeated events data.
#' 
#' @examples
#' 
#' # Load Golub & Steunenberg (2007) Data
#' data("GolubEUPData")
#' 
#' # Subset PURELY TO SPEED UP THE EXAMPLE
#' GolubEUPData <- GolubEUPData[1:1500, ]
#' 
#' # Expand data
#' GolubEUPDataExp <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
#'                      Time = 'begin', Time2 = 'end', event = 'event')
#' 
#' @seealso \code{\link{tvc}}
#' @importFrom data.table
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom DataCombine MoveFront
#' @importFrom DataCombine VarDrop
#' 
#' @export 

SurvExpand <- function(data, GroupVar, Time, Time2, event, messages = TRUE){
  # For CRAN
  FT <- UG <- FinalTime <- StartT <- ST <- FinishT <- NULL

  # Reminder
  if (isTRUE(messages)){
    message('\nReminder (1): currently SurvExpand does not support repeated events data.\n\n')
    message('Reminder (2): Times should be in 1 unit intervals or this could take awhile.\n\n')
  }
  # Warnings
  if (class(data[, Time]) != 'numeric'){
    if (isTRUE(messages)){
      message(paste0('Converting ', deparse(substitute(Time)), ' to numeric. Things might get wacky. Please check.'))
    }
    data[, Time] <- as.character(data[, Time])
    data[, Time] <- as.numeric(data[, Time])
  }
  if (class(data[, Time2]) != 'numeric'){
    if (isTRUE(messages)){
      message(paste0('Converting ', deparse(substitute(Time2)), ' to numeric. Things might get wacky. Please check.'))
    }
    data[, Time2] <- as.character(data[, Time2])
    data[, Time2] <- as.numeric(data[, Time2])
  }
  if (class(data[, event]) != 'numeric'){
    if (isTRUE(messages)){
      message(paste0('Converting ', deparse(substitute(event)), ' to numeric. Things might get wacky. Please check.'))
    }
    data[, event] <- as.character(data[, event])
    data[, event] <- as.numeric(data[, event])
  }

  # Create full unit-time data set
  Min <- min(data[, Time])
  Max <- max(data[, Time2])
  FullTimes = Min:Max
  FullTimes <- data.frame(FakeID = 1, FT = FullTimes)
  FullTimes <- data.table(FullTimes, key = "FakeID", allow.cartesian = TRUE)

  UniGroup <- unique(data[, GroupVar])
  UniGroup <- data.frame(FakeID = 1, UG = UniGroup)
  UniGroup <- data.table(UniGroup, key = "FakeID", allow.cartesian = TRUE)

  Full <- UniGroup[FullTimes, allow.cartesian = TRUE]
  Full <- Full[, !c('allow.cartesian.1', 'FakeID'), with = FALSE]
  Full <- setkey(Full, key = 'UG')

  if (isTRUE(messages)){
    message('\nExpanding data.\n') 
    message('Keeping only needed observations.\n')
  }

  # Keep only times needed for survival modelling
  End <- sort(unique(data[, Time2]))
  FullSub <- Full[FT %in% End]
  FullSub <- FullSub[order(UG, FT)]

  # Find last observation and subset data so that unit-observations greater than this are excluded
  DGroup <- eval(parse(text = paste0('group_by(data, ', GroupVar, ')')))
  DGroup <- eval(parse(text = paste0('dplyr::mutate(DGroup, FinalTime = max(', Time2, '))')))
  DGroup <- DGroup[, c(GroupVar, 'FinalTime')]
  names(DGroup) <- c('UG', 'FinalTime')
  DGroup <- DGroup[!duplicated(DGroup$UG, DGroup$FinalTime), ]
  DGroup <- data.table(DGroup, key = 'UG', allow.cartesian = TRUE)
  FullLast <- FullSub[DGroup, allow.cartesian = TRUE]
  FullLast <- FullLast[FT <= FinalTime]
  FullLast <- FullLast[, !c('allow.cartesian.1'), with = FALSE]

  # Create new time interval start variable
  FullLast <- FullLast[, 'ST' := FT - 1]

  # Merge with original data and create new time intervals
  DataMerge <- data
  Names <- c(GroupVar, Time, Time2, event)
  Names2 <- c('UG', 'StartT', 'FinishT', 'Event')
  DataMerge <- MoveFront(DataMerge, Names)
  colnames(DataMerge)[1:4] <- Names2 
  DataMerge <- data.table(DataMerge, key = 'UG', allow.cartisian = TRUE)

  FullComb <- FullLast[DataMerge, allow.cartesian = TRUE]
  FullComb <- FullComb[, !c('allow.cartisian', 'allow.cartesian'), with = FALSE]

  if (isTRUE(messages)){
    message('Doing a final clean up.')
  }

  Out <- FullComb[StartT <= ST & FinishT >= FT ]
  Out <- setkey(Out, NULL)
  Out <- unique(Out)

  # Final clean up
  Out <- data.frame(Out)
  Out$Event[Out$FinalTime != Out$FT] <- 0
  Out <- VarDrop(Out, c('StartT', 'FinishT', 'FinalTime'))
  Out <- MoveFront(Out, c('UG', 'ST', 'FT', 'Event'))
  colnames(Out)[1:4] <- Names

  return(Out)
}