
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qpNCA

<!-- a href='https://dplyr.tidyverse.org'><img src='man/figures/logo.png' align="right" height="139" /></a -->
<!-- badges: start 
[![Travis build status](https://travis-ci.com/qPharmetra/qpNCA.svg?branch=master)](https://travis-ci.com/qPharmetra/qpNCA)
[![Coverage status](https://codecov.io/gh/qpharmetra/qpNCA/branch/master/graph/badge.svg)](https://codecov.io/github/qpharmetra/qpNCA?branch=master)

 badges: end -->

## Overview

qpNCA is qPharmetra’s R package for Noncompartmental Analysis. It
provides regulatory-quality calculations of NCA parameters, using code
routinely tested against independently-validated output. You can find
appropriate manuals in the `vignettes` directory, and you can find tests
and reference data in the `tests` directory.

## Installation

The easiest way to get qpNCA is to install from CRAN:

``` r
install.packages("qpNCA")
```

Alternatively, install from github:

``` r
devtools::install_github('qpharmetra/qpNCA')
```

You can use the issues feature on
[github](https://github.com/qPharmetra/qpNCA/issues) to report issues.

## Usage

A good place to start is the help for `qpNCA()`.

``` r
library(qpNCA)
?qpNCA
example(qpNCA)
```

Or see the vignettes, like [this
one](http://htmlpreview.github.io/?https://github.com/qPharmetra/qpNCA/blob/master/vignettes/Example-stepwise-nca-analysis.html)
for example.
