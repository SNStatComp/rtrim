[![Build Status](https://travis-ci.org/markvanderloo/rtrim.svg)](https://travis-ci.org/markvanderloo/rtrim) 
[![Coverage Status](https://coveralls.io/repos/github/markvanderloo/rtrim/badge.svg?branch=master)](https://coveralls.io/github/markvanderloo/rtrim?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/rtrim)](http://cran.r-project.org/web/packages/rtrim)
[![Downloads](http://cranlogs.r-pkg.org/badges/rtrim)](http://cran.r-project.org/package=rtrim/) [![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)


# rtrim
Reimplementation of [TRIM](https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--) for R

### Installation and getting started

To install `rtrim`, issue the following command at the R console.
```r
install.packages("rtrim")
```

To get started, we recommend users to first work through the introductory vignette.This can be opened as follows.
```r
library(rtrim)
vignette("rtrim_for_TRIM_users")
```

An extended introduction showing all options can be opened as follows:
```r
library(rtrim)
vignette("Skylark_example")
```
After working through one of these vignettes, we advise users to have a look atthe help file of `rtrim`'s main function by typing `?trim`. Also, to get a feel of how the package works one can browse through the help files by following thelinks under the `see also` sections in each help file.

### Release planning

Version 1.0.2 will be uploaded to CRAN on 31 march 2017.

See the [NEWS](pkg/NEWS) file for a list of updates.


### Install the development version

Sometimes we release a beta version for public testing. 

**Current Beta Version:** 1.02


Beta versions can be installed as follows.

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
This procedure is necessary for the prerelease version. Once the package is on CRAN, installation of neither `Rtools` nor `drat` will be required.



### Note

The current package has been tested at two institutes and we have compared
numerical outcomes of hundreds of production jobs with outcomes of the original
TRIM software. We find that most results are reproducible to about `1E-4` or
better.  Nevertheless, since this software has not been tried and tested for
over 25 years like TRIM has, you can still expect to bump into a rare bug or
difference.

We would really appreciate it if you let us know using the following instructions.


### Telling us what works and what doesn't

We'd love to hear your findings. It would be great if you could do this by adding
issues to our [issues page](https://github.com/markvanderloo/rtrim/issues). You will 
need to create a github account for that. If you can not or do not want to create one,
send an e-mail to `rtrim AT cbs dot nl`.

When reporting bugs, it helps us a lot if you can create an  as small as possible example
that demonstrates the error. Some good general guidelines for reporting bugs are given [here](https://sifterapp.com/blog/2012/08/tips-for-effectively-reporting-bugs-and-issues/)






