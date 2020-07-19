test_that("coxsimtvc works for error in Issue # 24", {
    data("GolubEUPData")
    
    # Expand data (not run to speed processing time, but should be run)
    GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                               Time = 'begin', Time2 = 'end', event = 'event')
    
    # Create time interactions
    BaseVars <- c('qmv', 'qmvpostsea')
    GolubEUPData <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')
    
    # Run Cox PH Model
    M1 <- coxph(Surv(begin, end, event) ~ qmv + qmvpostsea,
                data = GolubEUPData, ties = "efron")
    
    # Create simtvc object for Relative Hazard 
    Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "qmv_log",
                      tfun = "log", from = 80, to = 2000,
                      Xj = 1, by = 15, ci = 0.99, nsim = 100)
    expect_gt(nrow(Sim1), 0) 
})
