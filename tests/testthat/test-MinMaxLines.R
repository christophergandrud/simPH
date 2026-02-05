test_that("MinMaxLines works with default byVars", {
    library(survival)

    data("CarpenterFdaData")

    # Run basic model
    M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
               deathrt1 + acutediz + hosp01  + hhosleng +
               mandiz01 + femdiz01 + peddiz01 + orphdum +
               vandavg3 + wpnoavg3 + condavg3 + orderent +
               stafcder, data = CarpenterFdaData)

    # Simulate Hazard Ratios
    Sim1 <- coxsimLinear(M1, b = "stafcder",
                         Xj = c(1237, 1600),
                         Xl = c(1000, 1000),
                         qi = "Hazard Ratio",
                         spin = TRUE, ci = 0.99,
                         nsim = 100)

    # Find summary statistics of the constricted interval
    Sum <- MinMaxLines(Sim1, clean = TRUE)

    # Test that the result is a data frame
    expect_true(is.data.frame(Sum))

    # Test that it has the expected columns
    expect_true("Xj" %in% names(Sum))
    expect_true("Min_CI" %in% names(Sum))
    expect_true("Lower50_CI" %in% names(Sum))
    expect_true("Median" %in% names(Sum))
    expect_true("Upper50_CI" %in% names(Sum))
    expect_true("Max_CI" %in% names(Sum))

    # Test that we get the expected number of rows (one per unique Xj value)
    expect_equal(nrow(Sum), 2)
})

test_that("MinMaxLines works with custom byVars", {
    library(survival)

    data("CarpenterFdaData")

    M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
               deathrt1 + acutediz + hosp01  + hhosleng +
               mandiz01 + femdiz01 + peddiz01 + orphdum +
               vandavg3 + wpnoavg3 + condavg3 + orderent +
               stafcder, data = CarpenterFdaData)

    Sim1 <- coxsimLinear(M1, b = "stafcder",
                         Xj = c(1237, 1600),
                         Xl = c(1000, 1000),
                         qi = "Hazard Ratio",
                         spin = TRUE, ci = 0.99,
                         nsim = 100)

    # Test with explicit byVars
    Sum <- MinMaxLines(Sim1, byVars = "Xj", clean = FALSE)

    expect_true(is.data.frame(Sum))
    expect_gt(nrow(Sum), 0)
})

test_that("MinMaxLines handles multiple grouping variables", {
    # Create test data with multiple grouping variables
    test_df <- data.frame(
        Xj = rep(c(1, 2), each = 100),
        Time = rep(c(10, 20), each = 50, times = 2),
        QI = rnorm(200, mean = 1, sd = 0.5)
    )

    # Test with multiple byVars
    Sum <- MinMaxLines(test_df, byVars = c("Xj", "Time"), clean = TRUE)

    expect_true(is.data.frame(Sum))
    # Should have 4 rows: 2 Xj values * 2 Time values
    expect_equal(nrow(Sum), 4)
    expect_true(all(c("Xj", "Time") %in% names(Sum)))
})
