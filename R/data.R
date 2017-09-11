# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Toy example -  observations for 5 correlated samples.
#'
#' A matrix containing the observations for the 1,000 positions of 5
#'  correlated samples.
#'@name toydata
#'@docType data
#' @format A simulated matrix with 1000 rows and 5 columns. Each column is a series
#' @keywords datasets
#' @examples
#' data(toyexample)
#' # Variational inference of a coupled hidden Markov Chains
#' resCHMM <- coupledHMM(X = toydata, nb.states = 3, S = cor(toystatus),
#'                       omega.list = c(0.3, 0.5, 0.7, 0.9))
#' # Breakpoints positions and status of segments
#' info <- clusterseg(resCHMM$status)
#' # head(info)
NULL

#' Toy example -  status for 5 correlated samples.
#'
#' A matrix containing the hidden status for the 1,000 positions of 5
#'  correlated samples.
#'@name toystatus
#'@docType data
#' @format A matrix of the hidden status corresponding to the \code{toydata} matrix.
#' @keywords datasets
NULL
