[![Build Status](https://travis-ci.org/markvanderloo/rtrim.svg)](https://travis-ci.org/markvanderloo/rtrim) 
[![Coverage Status](https://coveralls.io/repos/github/markvanderloo/rtrim/badge.svg?branch=master)](https://coveralls.io/github/markvanderloo/rtrim?branch=master)
# rtrim
Reimplementation of [TRIM](https://www.cbs.nl/en-gb/society/nature-and-environment/indices-and-trends--trim--) for R



### Install the beta version

A beta version can be installed from my [drat](https://cran.r-project.org/package=drat) repository. As this is a source distribution, you need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) if your OS is not unix-like (OSX, Linux, Solaris).
```
if(!require(drat)) install.packages('drat')
drat::addRepo('markvanderloo')
install.packages('rtrim',type='source')
```


