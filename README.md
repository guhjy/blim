# blim
Univariate Bayesian Linear Modelling made approachable for anyone: `blim()` works just like `lm()`!

## how to install

Has been tested on 64 bit R 3.2.5. Not guaranteed to work on other versions of R! Rtools is required for installation as a package. If you don't have it, you can download it from https://cran.r-project.org/bin/windows/Rtools/

```R
# install and load blim
library(devtools) # Required for installation
install_github("vankesteren/blim") # Install my package from github
library(blim) # Load this package
library(MASS) # forgot to add MASS in the package
# for usage:
?blim
```
