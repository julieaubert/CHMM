# CHMM

Coupled Hidden Markov Models


This package implements an exact and a variational inference for coupled Hidden Markov Models applied to the joint detection of copy number variations.


## Installation

You can install the released version of CHMM from [CRAN](https://CRAN.R-project.org) with:

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
## Reference

Wang, Xiaoqiang, Lebarbier, Emilie, Aubert, Julie and Robin, Stéphane. "Variational Inference for Coupled Hidden Markov Models Applied to the Joint Detection of Copy Number Variations" The International Journal of Biostatistics, vol. 15, no. 1, 2019, pp. 20180023. https://doi.org/10.1515/ijb-2018-0023

Aubert, Julie, Wang, Xiaoqiang, Lebarbier, Emilie, Robin, Stéphane. "CHMM : an R package for coupled Hidden Markov Models". R User Conference 2017, Jul 2017, Bruxelles, Belgium. ⟨hal-01661257⟩

