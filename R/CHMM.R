# CHMM R package
# Copyright INRA 2017
# UMR MIA-Paris, AgroParisTech, INRA, Universite Paris-Saclay, 75005, Paris, France
###################################################################
#' Hidden status of 10 correlated samples.
#'
#' A matrix containing the hidden status fo the 1,000 positions of 10
#'  correlated samples.
#'@name toyexample
#'@docType data
#' @format A list containing a matrix with 1000 rows and 10 columns:
#' \describe{
#'   \item{X}{a matrix with 1000 rows and 5 columns. Each column is a series}
#'   \item{hidden}{a matrix of the corresponding hidden status}
#' }
#' @keywords datasets
#' @export
#' @examples
#' data(toyexample)
#' # Variational inference of a coupled hidden Markov Chains
#' resCHMM <- chmm(X, nb.states = 3, S = cor(hidden), omega.list = c(0.3, 0.5, 0.7, 0.9))

#' A matrix containing the hidden status fo the 1,000 positions of 10
#'  correlated samples.
#'@name toyexample
#'@docType data
#' @format A list containing a matrix with 1000 rows and 10 columns:
#' \describe{
#'   \item{X}{a matrix with 1000 rows and 5 columns. Each column is a series}
#'   \item{hidden}{a matrix of the corresponding hidden status}
#' }
#' @keywords datasets
#' @export
#' @examples
#' data(toyexample)
#' # Variational inference of a coupled hidden Markov Chains
#' resCHMM <- chmm(X, nb.states = 3, S = cor(hidden), omega.list = c(0.3, 0.5, 0.7, 0.9))
