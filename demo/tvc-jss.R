#############
# Time-interactive effect examples from JSS paper
# Christopher Gandrud
# Updated 5 April 2014
#############

# Load packages
library("survival")
library("simPH")
library("ggplot2")
library("gridExtra")


##### Illustration of time-varying interactive effects ######
# Load Golub & Steunenberg (2007) data. The data is included with simPH.
data("GolubEUPData")

# Examine the data's format
head(GolubEUPData[, 2:5])

# Expand data into equally spaced time intervals
 GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                      Time = 'begin', Time2 = 'end', event = 'event') 

# Examine the data's new format
head(GolubEUPData[, 1:4])

# Create natural log-time interactions
Golubtvc <- function(x){
  tvc(data = GolubEUPData, b = x, tvar = "end", tfun = "log")
}
GolubEUPData$Lqmv <- Golubtvc("qmv")
GolubEUPData$Lbacklog <- Golubtvc("backlog")
GolubEUPData$Lcoop <- Golubtvc("coop")
GolubEUPData$Lcodec <- Golubtvc("codec")
GolubEUPData$Lqmvpostsea <- Golubtvc("qmvpostsea")
GolubEUPData$Lthatcher <- Golubtvc("thatcher")

# Estimate model
M2 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea + qmvpostteu +
              coop + codec + eu9 + eu10 + eu12 + eu15 + thatcher + 
              agenda + backlog + Lqmv + Lqmvpostsea + Lcoop + Lcodec +
              Lthatcher + Lbacklog,
            data = GolubEUPData, ties = "efron")

## Create simtvc object for first difference (central interval)
Sim3_1 <- coxsimtvc(obj = M2, b = "qmv", btvc = "Lqmv",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000,
                    by = 5, ci = 0.95)

# Create simtvc object for first difference (SPIn)
Sim3_2 <- coxsimtvc(obj = M2, b = "qmv", btvc = "Lqmv",
                    qi = "First Difference", Xj = 1,
                    tfun = "log", from = 80, to = 2000,
                    by = 5, ci = 0.95, spin = TRUE)

# Create first difference plots
Plot3_1 <- simGG(Sim3_1, xlab = "\nTime in Days", 
                 title = "Central Interval\n", alpha = 0.3,
                 type = "ribbons", lsize = 0.5, legend = FALSE)

Plot3_2 <- simGG(Sim3_2, ylab = "", xlab = "\nTime in Days",
                 title = "SPIn\n", alpha = 0.3,
                 type = "ribbons", lsize = 0.5, legend = FALSE)

# Combine plots
grid.arrange(Plot3_1, Plot3_2, ncol = 2)

# Create simtvc object for relative hazard
Sim4 <- coxsimtvc(obj = M2, b = "backlog", btvc = "Lbacklog",
                  qi = "Relative Hazard", Xj = seq(40, 200, 40),
                  tfun = "log", from = 1200, to = 7000, by = 100,
                  nsim = 200)

# Create relative hazard plot
simGG(Sim4, xlab = "\nTime in Days", type = "ribbons",
      leg.name = "Backlogged \n Items")
