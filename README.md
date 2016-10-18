[![Build Status](https://travis-ci.org/markvanderloo/rtrim.svg)](https://travis-ci.org/markvanderloo/rtrim) 
[![Coverage Status](https://coveralls.io/repos/github/markvanderloo/rtrim/badge.svg?branch=master)](https://coveralls.io/github/markvanderloo/rtrim?branch=master)

# rtrim
Reimplementation of [TRIM](https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--) for R



### Install the review version (beta)


A review release (version 0.0.3.1) can now be installed with the following instructions.


1. If you are a Windows user, first install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
2. Open an R session, and issue the following commands.
```r
# we need the 'drat' package to install the review version:
if(!require(drat)) install.packages('drat')
# we tell R where to find the review version (on github)
drat::addRepo('markvanderloo')
# install as usual
install.packages('rtrim',type='source')
```
This procedure is necessary while we are still in beta. Once the package is on CRAN, installation neither `Rtools` nor `drat` will be required.

#### What is in the review version?

The current version has most, but not all functionality of the original TRIM
software. In short, the following features are currently supported.

- Estimate coefficients for TRIM model 2
- Estimate coefficients TRIM model 3
- Use of covariates
- Estimating serial correlation 
- Estimating overdispersion
- Manual changepoint specification (model 2)
- Stepwise changepoint selection (model 2)
- Goodness-of-fit: AIC, p-value and Xi-squared
- Wald test statistic for significance of slopes
- Computing time totals
- Computing indices
- Computing overall slope (slope over multiple pieces of the piecewise linear model)
- read instructions TRIM command files (for backward compatibility with TRIM)
- read TRIM data files (for backward compatibility with TRIM)


#### What is not (yet) in the review version?

- Estimate coefficients for TRIM model 1
- The `AUTODELETE` function has known problems
- Writing, reading and reusing covariance matrices


### Getting started

To get started, we recommend users to first work through the introductory vignette.
This can be opened as follows.
```r
library(rtrim)
vignette("rtrim_for_TRIM_users")
```

After working through the vignette, we advise users to have a look at the help
file of `rtrim`'s main function by typing `?trim`. Also, to get a feel of how
the package works one can browse through the help files by following the links
under the `see also` sections in each help file.




### Telling us what works and what doesn't

We'd love to hear your findings. It would be great if you could do this by adding
issues to our [issues page](https://github.com/markvanderloo/rtrim/issues). You will 
need to create a github account for that. If you can not or do not want to create one,
send an e-mail to `mark dot vanderloo @ gmail dot com`.

When reporting bugs, it helps us a lot if you can create an  as small as possible example
that demonstrates the error. Some good general guidelines for reporting bugs are given [here](https://sifterapp.com/blog/2012/08/tips-for-effectively-reporting-bugs-and-issues/)










