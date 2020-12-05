test_that("coxsimtvc works for error in Issue # 24", {
    
    # Hack to address https://github.com/christophergandrud/simPH/issues/26
    # Modified from https://conjugateprior.org/2015/06/identifying-the-os-from-r/
    # get_os <- function(){
    #     sysinf <- Sys.info()
    #     if (!is.null(sysinf)){
    #         os <- sysinf['sysname']
    #         if (os == 'Darwin')
    #             os <- "macos"
    #     } else { ## mystery machine
    #         os <- .Platform$OS.type
    #         if (grepl("^darwin", R.version$os))
    #             os <- "macos"
    #         if (grepl("linux-gnu", R.version$os))
    #             os <- "linux"
    #     }
    #     tolower(os)
    # }
    
    
    data("GolubEUPData")
    
    # Expand data (not run to speed processing time, but should be run)
#    if (get_os() != "linux") {
        GolubEUPData <- SurvExpand(GolubEUPData, GroupVar = 'caseno',
                                   Time = 'begin', Time2 = 'end', event = 'event')
        
        # Create time interactions
        BaseVars <- c('qmv')
        GolubEUPData <- tvc(GolubEUPData, b = BaseVars, tvar = 'end', tfun = 'log')
        
        # Run Cox PH Model
        M1 <- coxph(Surv(begin, end, event) ~ qmv + qmv_log,
                    data = GolubEUPData, ties = "efron")
        
        # Create simtvc object for Relative Hazard 
        Sim1 <- coxsimtvc(obj = M1, b = "qmv", btvc = "qmv_log",
                          tfun = "log", from = 80, to = 2000,
                          Xj = 1, by = 15, ci = 0.99, nsim = 100)
        expect_gt(nrow(Sim1), 0) 
#    }
})
