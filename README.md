<img src="vignettes/logo-epiphy-01.png" alt="epiphy logo"/>

# epiphy: An R package to analyze plant disease epidemics

**epiphy** provides a common framework for spatialized plant disease intensity data collected at one or more time points. Several statistical methods to describe and quantify plant disease epidemics are implemented in this package.

## Resources

The online documentation is available at [epiphy.org](http://epiphy.org).

## Installation

Install the latest development version from **GitHub**:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
    # ^ Returns FALSE if devtools is not installed.
    install.packages("devtools")
}
devtools::install_github("chgigot/epiphy")
```

**Note:** You can also use the following command to build the vignettes when installing the package:

```r
devtools::install_github("chgigot/epiphy", build_vignettes = TRUE)
```
