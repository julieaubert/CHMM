---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# CHMM

Coupled Hidden Markov Models


This package implements an exact and a variational inference for coupled Hidden Markov Models applied to the joint detection of copy number variations.

## Installation

You can install the released version of cobiclust from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CHMM")
```

## Example

This is a toy example (1 000 positions of 5 correlated samples) which shows you how to use the package.

```{r example}
data(toyexample)
# Variational inference of a coupled hidden Markov Chains
resCHMM <- coupledHMM(X = toydata, nb.states = 3, S = cor(toystatus),
                      omega.list = c(0.3, 0.5, 0.7, 0.9))
# Breakpoints positions and status of segments
info <- clusterseg(resCHMM$status)
```

