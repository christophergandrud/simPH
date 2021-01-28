# simPH

<img src="man/figures/plotbanner.png" height="40" width="1000" alt="banner-image"></img>

**Christopher Gandrud**

### [![R build status](https://github.com/christophergandrud/simPH/workflows/R-CMD-check/badge.svg)](https://github.com/christophergandrud/simPH/actions) [![Build Status](https://travis-ci.org/christophergandrud/simPH.png)](https://travis-ci.org/christophergandrud/simPH) [![Codecov test coverage](https://codecov.io/gh/christophergandrud/simPH/branch/master/graph/badge.svg)](https://codecov.io/gh/christophergandrud/simPH?branch=master) [![CRAN Version](https://CRAN.R-project.org/package=simPH)](https://CRAN.R-project.org/package=simPH) ![CRAN Downloads](https://cranlogs.r-pkg.org/badges/last-month/simPH) ![CRAN Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/simPH)

#### Please report any bugs to:

<https://github.com/christophergandrud/simPH/issues>

---

**simPH** is an R package for simulating and plotting quantities of interest
(relative hazards, first differences, and hazard ratios) for linear
coefficients, multiplicative interactions, polynomials, penalised splines, and
non-proportional hazards, as well as stratified survival curves from Cox
Proportional Hazard models.

For more information plus examples, please see the
[description paper](https://www.jstatsoft.org/article/view/v065i03) in the
[Journal of Statistical Software](https://www.jstatsoft.org/).

To cite the paper please use:

```
@article{simPH_JSS,
    author = {Christopher Gandrud},
    title = {simPH: An R Package for Illustrating Estimates from Cox
        Proportional Hazard Models Including for Interactive and Nonlinear
        Effects},
    journal = {Journal of Statistical Software},
    year = {2015},
    volume = {65},
    issue = {3},
    pages = {1--20}
}
```

## Functions

The package includes the following functions:

#### Simulation Functions

- `coxsimLinear`: a function for simulating relative hazards, first differences,
hazard ratios, and hazard rates for linear, non-time interacted covariates from
Cox Proportional Hazard models.

- `coxsimtvc`: a function for simulating time interactive hazards (relative
hazards, first differences, and hazard ratios) for covariates from Cox
Proportional Hazard models. The function will calculate time-interactive hazard
ratios for multiple strata estimated from a stratified Cox Proportional Hazard
model.

- `coxsimSpline`: a function for simulating quantities of interest from
penalised splines using multivariate normal distributions. Currently does not
support simulating hazard rates from stratified models. **Note:** be extremely
careful about the number of simulations you ask the function to find. It is very
easy to ask for more than your computer can handle.

- `coxsimPoly`: a function for simulating quantities of interest for a range of
values for a polynomial nonlinear effect from Cox Proportional Hazard models.

- `coxsimInteract`: a function for simulating quantities of interest for linear
multiplicative interactions, including marginal effects and hazard rates.

#### Plotting Functions

Results from these functions can be plotted using the `simGG` method. The
syntax and capabilities of `simGG` varies depending on the `sim` object class
you are using:

- `simGG.simlinear`: plots simulated linear time-constant hazards using
[ggplot2](https://CRAN.R-project.org/package=ggplot2 ).

- `simGG.simtvc`: uses **ggplot2** to graph the simulated time-varying relative
hazards, first differences, hazard ratios or stratified hazard rates.

- `simGG.simspline`: uses **ggplot2** to plot
quantities of interest from `simspline` objects, including relative hazards,
first differences, hazard ratios, and hazard rates.

- `simGG.simpoly`: uses **ggplot2** to graph the simulated polynomial quantities
of interest.

- `simGG.siminteract`: uses **ggplot2** to graph linear multiplicative
interactions.

##### Additional styling

Because in almost all cases `simGG` returns a *ggplot2* object, you can add
additional aesthetic attributes in the normal *ggplot2* way. See the
[ggplot2 documentation for more details](https://ggplot2.tidyverse.org/reference/).

#### Misc.

- `SurvExpand`: a function for converting a data frame of non-equal interval
continuous observations into equal interval continuous observations. This is
useful to do before creating time interactions.

    + See also the [SurvSetup](https://github.com/christophergandrud/SurvSetup)
    package for additional functions to set up your data for survival analysis.

- `tvc`: a function for creating time interactions. Currently supports
`'linear'`, natural `'log'`, and exponentiation (`'power'`).

- `setXl`: a function for setting valid `Xl` values given a sequence of fitted
`Xj` values. This makes it more intuitive to find hazard ratios and first
differences for comparisons between some Xj fitted values and Xl values other
than 0.

- `ggfitStrata`: a function to plot fitted stratified survival curves estimated
from `survfit` using **ggplot2**. This function builds on the **survival**
package's `plot.survfit` function. One major advantage is the ability to split
the survival curves into multiple plots and arrange them in a grid. This makes
it easier to examine many strata at once. Otherwise they can be very bunched up.

- `MinMaxLines`: a function for summarising the constricted intervals from the
simulations, including the median, upper and lower bounds and
the middle 50% of these intervals.

## Installation

The package is available on CRAN and can be installed in the normal R way.

## Tip

Before running the simulation and graph functions in this package carefully
consider how many simulations you are about to make. Especially for hazard rates
over long periods of time and with multiple strata, you can be asking **simPH**
to run very many simulations. This will be computationally intensive.

## Sources

### Simulating Parameter Estimates

For more information about simulating parameter estimates to make interpretation
of results easier see:

Licht, Amanda A. 2011. [“Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.”](https://www.jstor.org/stable/23011265)
*Political Analysis* 19: 227–43.

King, Gary, Michael Tomz, and Jason Wittenberg. 2000.
[“Making the Most of Statistical Analyses: Improving Interpretation and Presentation.”](https://www.jstor.org/stable/2669316) *American Journal of Political Science* 44(2): 347–61.

### Stratified Cox PH

For more information about stratified Cox PH models see:

Box-Steffensmeier, Janet M, and Suzanna De Boef. 2006. [“Repeated Events Survival Models: the Conditional Frailty Model.”](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2434) *Statistics in Medicine* 25(20): 3518–33.

### Shortest Probability Intervals

To learn more about shortest probability intervals (and also for the source of
the code that made this possible in **simPH**) see:

Liu, Y., Gelman, A., & Zheng, T. (2015).
["Simulation-efficient Shortest Probablility Intervals."](https://www.stat.columbia.edu/~gelman/research/published/spin.pdf)
*Statistics and Computing* 25:809-819.

**Also good:** Hyndman, R. J. (1996).
["Computing and Graphing Highest Density Regions."](https://www.jstor.org/stable/10.2307/2684423)
*The American Statistician*, 50(2): 120–126.

### Interpreting Interactions

For more information about interpreting interaction terms:

Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006.
[“Understanding Interaction Models: Improving Empirical Analyses.”](https://www.jstor.org/stable/25791835)
*Political Analysis* 14(1): 63–82.

### The Olden Days

For an example of how non-proportional hazard results were often presented
before **simPH** see (some of the problems I encountered in this paper were a
major part of why I developed this package):

Gandrud, Christopher. 2013. [“The Diffusion of Financial Supervisory Governance Ideas.”](https://www.tandfonline.com/doi/full/10.1080/09692290.2012.727362)
*Review of International Political Economy*. 20(4): 881-916.

## Thoughts on future versions

I intended to expand the quantities of interest that can be simulated and graphed
for Cox PH models. I am also currently working on functions that can simulate
and graph hazard ratios estimated from
[Fine and Gray competing risks models](https://www.jstor.org/stable/2670170).

It would also be great to enable graphing hazard ratios with frailties.

However, I don't currently plan to make major improvements to the package. My ambition is to ensure that it remains stable and on CRAN.

But pull requests are always welcome!

---

Licensed under [GPL-3](https://github.com/christophergandrud/simPH/blob/master/LICENSE.md)
