# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
## __________________________________________________________
##
## Function :: Emis.Gauss()
## __________________________________________________________
##
#' Emis.Gauss
#'
#' @param X data matrix of observations.
#' @param esAvg a numeric of the estimated mean for each state.
#' @param esVar a numeric of the estimated variance for each state.
#' @importFrom stats dnorm
#' @keywords internal
#'
Emis.Gauss <- function(X, esAvg, esVar){
  nbS <- length(esAvg)
  apply(as.matrix(1:nbS), 1, function(r) dnorm(X, mean = esAvg[r], sd = sqrt(esVar[r])) )
}
## __________________________________________________________
##
## Function :: Emis.LRR()
## __________________________________________________________
##

#' Emis.LRR
#'
#' @param X data matrix of observations.
#' @param esAvg a numeric of the estimated mean for each state.
#' @param esVar a numeric of the estimated variance for each state.
#' @param weight weight.
#' @importFrom stats dnorm
#' @keywords internal
#'
#'
Emis.LRR <- function(X, esAvg, esVar, weight = 0.05){
  nbS <- length(esAvg)
  apply(as.matrix(1:nbS), 1, function(r)
    weight + (1 - weight) * dnorm(X, mean = esAvg[r], sd = sqrt(esVar[r])) )
}
## __________________________________________________________
##
## Function :: EmisGb.Gauss()
## __________________________________________________________
##
#' EmisGb.Gauss
#'
#' @param X data matrix of observations.
#' @param esAvgGb a matrix of \code{nbK} rows and \code{nbI} columns of estimated mean.
#' @param esVarGb a matrix of \code{nbK} rows and \code{nbI} columns of estimated variance.
#' @importFrom stats dnorm
#' @keywords internal
#'
EmisGb.Gauss <- function(X, esAvgGb, esVarGb){
  nbK <- nrow(esAvgGb)
  emisGb <- matrix(0, nrow(X), nbK)
  for(k in 1:nbK) {
    emisTmp <- NULL
    for(i in 1:ncol(X)) {
      emisTmp <- cbind(emisTmp,  dnorm(X[,i], mean = esAvgGb[k,i], sd = sqrt(esVarGb[k,i])) )
    }
    emisGb[,k] <- exp(rowSums(log(emisTmp)))
  }
  emisGb
}

