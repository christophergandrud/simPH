simPH
======

### Christopher Gandrud

### Version 0.8.4

### Please report any bugs at <https://github.com/christophergandrud/simPH/issues>.

---

An R package for simulating and plotting quantities of interest (relative hazards, first differences, and hazard ratios) for linear coefficients, multiplicative interactions, polynomials, penalised splines, and non-proportional hazards, as well as stratified survival curves from Cox Proportional Hazard models.

For more information plus examples, please see this [working paper](http://ssrn.com/abstract=2318977).

The package includes the following functions:

#### Simulation Functions

- `coxsimLinear`: Simulates relative hazards, first differences, hazard ratios, and hazard rates for linear time-constant covariates from Cox Proportional Hazard models.

- `coxsimtvc`: a function for simulating time-varying hazards (relative hazards, first differences, and hazard ratios) from a Cox PH model estimated using `coxph` from the [survival](http://cran.r-project.org/web/packages/survival/index.html) package. For more information see this [blog post](http://christophergandrud.blogspot.kr/2012/10/graphing-non-proportional-hazards-in-r.html). The function will calculate time-varying hazard ratios for multiple strata estimated from a stratified Cox PH model.

- `coxsimSpline`: a function for simulating quantities of interest from penalised splines using multivariate normal distributions. Currently does not support simulating hazard rates from stratified models. **Note:** be extremely careful about the number of simulations you ask the function to find. It is very easy to ask for more than your computer can handle.

- `coxsimPoly`: a function for simulating quantities of interest for a range of values for a polynomial nonlinear effect from Cox Proportional Hazards models.

- `coxsimInteract`: a function for simulating quantities of interest for linear multiplicative interactions, including marginal effects and hazard rates.

#### Plotting Functions

**Important Change:** From **simPH** version 0.05 you can now use the method `simGG` for plotting all simulations rather than the old sim model object specific commands. In practical terms this means that you can now just use the command `simGG` rather than the old gg. . . commands. The old commands are deprecated and will no longer work.

The syntax and capabilities of `simGG` varies depending on the sim object class you are using:

- `simGG.simlinear`: plots simulated linear time-constant hazards using [ggplot2](http://ggplot2.org/).

- `simGG.simtvc`: uses **ggplot2** to graph the simulated time-varying relative hazards, first differences, hazard ratios or stratified hazard rates.

- `simGG.simspline`: uses **ggplot2** and `scatter3d` (from the [car](http://cran.r-project.org/web/packages/car/index.html) package) to plot quantities of interest from `simspline` objects, including relative hazards, first differences, hazard ratios, and hazard rates.

- `simGG.simpoly`: uses **ggplot2** to graph the simulated polynomial quantities of interest.

- `simGG.siminteract`: uses **ggplot2** to graph linear multiplicative interactions.

##### Additional styling

Because in almost all cases `simGG` returns a *ggplot2* object, you can add additional aesthetic attributes in the normal *ggplot2* way. See the [ggplot2 documentation for more details](http://docs.ggplot2.org/current/). Here is a quick example:


---

#### Misc.

- `tvc`: a function for creating time interactions. Currently supports `'linear'`, natural `'log'`, and exponentiation (`'power'`).

- `ggfitStrata`: a function to plot fitted stratified survival curves estimated from `survfit` using **ggplot2**. This function builds on the **survival** package's `plot.survfit` command. One major advantage is the ability to split the survival curves into multiple plots and arrange them in a grid. This makes it easier to examine many strata at once. Otherwise they can be very bunched up.

## Installation

The package is available on CRAN and can be installed in the normal R way.

To install the development version use the [devtools](https://github.com/hadley/devtools) command `install_github`. Here is the code for installing the most recent development version:

```r
devtools::install_github("simPH", "christophergandrud")
```

## Tip

Before running the simulation and graph commands in this package carefully consider how many simulations you are about to make. Especially for hazard rates over long periods of time and with multiple strata, you can be asking **simPH** to run very many simulations. This will be computationally intensive. 

---

## Sources

### Simulating Parameter Estimates

For more information about simulating parameter estimates to make interpretation of results easier see:

Licht, Amanda A. 2011. [“Change Comes with Time: Substantive Interpretation of Nonproportional Hazards in Event History Analysis.”](http://pan.oxfordjournals.org/content/19/2/227.abstract) Political Analysis 19: 227–43.

King, Gary, Michael Tomz, and Jason Wittenberg. 2000. [“Making the Most of Statistical Analyses: Improving Interpretation and Presentation.”](http://www.jstor.org/stable/2669316) American Journal of Political Science 44(2): 347–61.

### Stratified Cox PH

For more information about stratified Cox PH models (and frailties, which I am working to incorporate in future versions) see:

Box-Steffensmeier, Janet M, and Suzanna De Boef. 2006. [“Repeated Events Survival Models: the Conditional Frailty Model.”](http://onlinelibrary.wiley.com/doi/10.1002/sim.2434/abstract;jsessionid=28218243DD3D6E01A3D10EEE75D96675.d01t02) Statistics in Medicine 25(20): 3518–33.

### Shortest Probability Intervals

To learn more about shortest probability intervals (and also for the source of the code that made this possible in **simPH**) see:

Liu, Y., Gelman, A., & Zheng, T. (2013). ["Simulation-efficient Shortest Probablility Intervals."](http://arxiv.org/pdf/1302.2142v1.pdf) Arvix. 

**Also good:** Hyndman, R. J. (1996). ["Computing and Graphing Highest Density Regions."](http://www.jstor.org/stable/10.2307/2684423) The American Statistician, 50(2), 120–126.

### Interpreting Interactions

For more information about interpreting interaction terms:

Brambor, Thomas, William Roberts Clark, and Matt Golder. 2006. [“Understanding Interaction Models: Improving Empirical Analyses.”](http://pan.oxfordjournals.org/content/14/1/63.abstract) Political Analysis 14(1): 63–82.

### The Olden Days

For an example of how non-proportional hazard results were often presented before **simPH** see (some of the problems I encountered in this paper were a major part of why I'm developing this package): 

Gandrud, Christopher. 2013. [“The Diffusion of Financial Supervisory Governance Ideas.”](http://www.tandfonline.com/doi/full/10.1080/09692290.2012.727362) Review of International Political Economy. 20(4):881-916.

---

## Future Plans

I intend to expand the quantities of interest that can be simulated and graphed for Cox PH models. I am also currently working on functions that can simulate and graph hazard ratios estimated from [Fine and Gray competing risks models](http://www.jstor.org/stable/2670170). 

I am also working on a way to graph hazard ratios with frailties. 
