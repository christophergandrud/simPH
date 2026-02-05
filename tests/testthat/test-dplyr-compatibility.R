# Tests for dplyr compatibility after removing deprecated functions
# These tests verify that the migration from group_by_(), distinct_(),
# group_by_at(), and as_data_frame() to modern dplyr functions works correctly.

test_that("coxsimLinear works with modern dplyr (tests IntervalConstrict)", {
    library(survival)

    data("CarpenterFdaData")

    M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
               deathrt1 + acutediz + hosp01  + hhosleng +
               mandiz01 + femdiz01 + peddiz01 + orphdum +
               vandavg3 + wpnoavg3 + condavg3 + orderent +
               stafcder, data = CarpenterFdaData)

    # Test Relative Hazard (uses IntervalConstrict internally)
    Sim1 <- coxsimLinear(M1, b = "stafcder",
                         Xj = c(1237, 1600),
                         Xl = c(1000, 1000),
                         qi = "Relative Hazard",
                         ci = 0.95,
                         nsim = 100)

    # coxsim objects are lists with a $sims data.frame component
    expect_s3_class(Sim1, "simlinear")
    expect_s3_class(Sim1, "coxsim")
    expect_true(is.data.frame(as.data.frame(Sim1)))
    expect_gt(nrow(as.data.frame(Sim1)), 0)
})

test_that("coxsimLinear with Hazard Ratio works (tests IntervalConstrict with spin)", {
    library(survival)

    data("CarpenterFdaData")

    M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
               deathrt1 + acutediz + hosp01  + hhosleng +
               mandiz01 + femdiz01 + peddiz01 + orphdum +
               vandavg3 + wpnoavg3 + condavg3 + orderent +
               stafcder, data = CarpenterFdaData)

    # Test with spin = TRUE which uses different IntervalConstrict path
    Sim1 <- coxsimLinear(M1, b = "stafcder",
                         Xj = c(1237, 1600),
                         Xl = c(1000, 1000),
                         qi = "Hazard Ratio",
                         spin = TRUE,
                         ci = 0.99,
                         nsim = 100)

    expect_s3_class(Sim1, "coxsim")
    expect_gt(nrow(as.data.frame(Sim1)), 0)

    # This is the exact test case from the CRAN error
    Sum <- MinMaxLines(Sim1, clean = TRUE)
    expect_true(is.data.frame(Sum))
    expect_equal(nrow(Sum), 2)
})

test_that("coxsimLinear with First Difference works", {
    library(survival)

    data("CarpenterFdaData")

    M1 <- coxph(Surv(acttime, censor) ~ prevgenx + lethal +
               deathrt1 + acutediz + hosp01  + hhosleng +
               mandiz01 + femdiz01 + peddiz01 + orphdum +
               vandavg3 + wpnoavg3 + condavg3 + orderent +
               stafcder, data = CarpenterFdaData)

    # Test First Difference qi
    Sim1 <- coxsimLinear(M1, b = "stafcder",
                         Xj = c(1237, 1600),
                         Xl = c(1000, 1000),
                         qi = "First Difference",
                         ci = 0.95,
                         nsim = 100)

    expect_s3_class(Sim1, "coxsim")
    expect_gt(nrow(as.data.frame(Sim1)), 0)
})
